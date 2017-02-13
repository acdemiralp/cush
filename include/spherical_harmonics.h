#ifndef CUSH_SPHERICAL_HARMONICS_H_
#define CUSH_SPHERICAL_HARMONICS_H_

#define _USE_MATH_DEFINES

#include <math.h>

#include <device_launch_parameters.h>
#include <vector_types.h>

#include <clebsch_gordan.h>
#include <decorators.h>
#include <legendre.h>

// Based on "Spherical Harmonic Lighting: The Gritty Details" by Robin Green.
namespace cush
{
INLINE COMMON unsigned int maximum_degree   (const unsigned int coefficient_count)
{
  return sqrtf(coefficient_count) - 1;
}
INLINE COMMON unsigned int coefficient_count(const unsigned int max_l)
{
  return powf(max_l + 1, 2);
}
INLINE COMMON unsigned int coefficient_index(const unsigned int l, const int m)
{
  return l * (l + 1) + m;
}
INLINE COMMON int2         coefficient_lm   (const unsigned int index)
{
  int2 lm;
  lm.x = floor(sqrtf(index));
  lm.y = index - powf(lm.x, 2) - lm.x;
  return lm;
}

template<typename precision>
COMMON precision evaluate(
  const unsigned int l    ,
  const          int m    ,
  const precision&   theta,
  const precision&   phi  )
{
  precision kml = sqrt((2.0 * l + 1) * factorial<precision>(l - abs(m)) / 
                       (4.0 * M_PI   * factorial<precision>(l + abs(m))));
  if (m > 0)
    return sqrt(2.0) * kml * cos( m * theta) * associated_legendre(l,  m, cos(phi));
  if (m < 0)
    return sqrt(2.0) * kml * sin(-m * theta) * associated_legendre(l, -m, cos(phi));
  return kml * associated_legendre(l, 0, cos(phi));
}
template<typename precision>
COMMON precision evaluate(
  const unsigned int index,
  const precision&   theta,
  const precision&   phi  )
{
  auto lm = coefficient_lm(index);
  return evaluate(lm.x, lm.y, theta, phi);
}

// Not used internally as the two for loops also need to be parallelized.
template<typename precision>
COMMON precision evaluate_sum(
  const unsigned int max_l       ,
  const precision*   coefficients,
  const precision&   theta       ,
  const precision&   phi         )
{
  precision sum = 0.0;
  for (int l = 0; l <= max_l; l++)
    for (int m = -l; m <= l; m++)
      sum += evaluate(l, m, theta, phi) * coefficients[coefficient_index(l, m)];
  return sum;
}

// Based on "Rotation Invariant Spherical Harmonic Representation of 3D Shape Descriptors" by Kazhdan et al.
template<typename precision>
COMMON precision compare(
  const unsigned int coefficient_count,
  const precision*   lhs_coefficients ,
  const precision*   rhs_coefficients )
{
  precision value = 0;
  for (auto index = 0; index < coefficient_count; index++)
    value += pow(lhs_coefficients[index] - rhs_coefficients[index], 2);
  return sqrt(value);
}

// Call on a vectors_size x coefficient_count(max_l) 2D grid.
template<typename vector_type, typename precision>
GLOBAL void calculate_matrix(
  const unsigned int max_l        ,
  const unsigned int vectors_size ,
  const vector_type* vectors      , // 1D: vectors_size
  precision*         output_matrix) // 2D: vectors_size * coefficient_count(max_l)
{
  auto vector_index      = blockIdx.x * blockDim.x + threadIdx.x;
  auto coefficient_index = blockIdx.y * blockDim.y + threadIdx.y;
  
  if (vector_index      > vectors_size || 
      coefficient_index > coefficient_count(max_l))
    return;

  auto& vector = vectors[vector_index];
  
  // Note: Column first indexing.
  atomicAdd(&output_matrix[vector_index + vectors_size * coefficient_index], evaluate(coefficient_index, vector.y, vector.z));
}
// Call on a dimensions.x x dimensions.y x dimensions.z 3D grid.
template<typename vector_type, typename precision>
GLOBAL void calculate_matrices(
  const uint3        dimensions     ,
  const unsigned int max_l          ,
  const unsigned int vectors_size   ,
  const vector_type* vectors        , // 4D: dimensions.x * dimensions.y * dimensions.z * coefficient_count(max_l)
  precision*         output_matrices) // 5D: dimensions.x * dimensions.y * dimensions.z * vectors_size * coefficient_count(max_l)
{
  auto x = blockIdx.x * blockDim.x + threadIdx.x;
  auto y = blockIdx.y * blockDim.y + threadIdx.y;
  auto z = blockIdx.z * blockDim.z + threadIdx.z;
  
  if (x > dimensions.x || y > dimensions.y || z > dimensions.z)
    return;
  
  auto volume_index    = z + dimensions.z * (y + dimensions.y * x);
  auto vectors_offset  = volume_index * vectors_size;
  auto matrices_offset = volume_index * vectors_size * coefficient_count(max_l);

  calculate_matrix<<<dim3(vectors_size, coefficient_count(max_l)), 1>>>(
    max_l, 
    vectors_size, 
    vectors         + vectors_offset , 
    output_matrices + matrices_offset);
}

// Call on a output_resolution.x x output_resolution.y x coefficient_count(max_l) 3D grid.
template<typename precision, typename point_type>
GLOBAL void sample_sum(
  const unsigned int max_l            ,
  const uint2        output_resolution,
  const precision*   coefficients     ,
  point_type*        output_points    ,
  unsigned int*      output_indices   )
{
  auto longitude_index   = blockIdx.x * blockDim.x + threadIdx.x;
  auto latitude_index    = blockIdx.y * blockDim.y + threadIdx.y;
  auto coefficient_index = blockIdx.z * blockDim.z + threadIdx.z;
  
  if (longitude_index   > output_resolution.x    ||
      latitude_index    > output_resolution.y    ||
      coefficient_index > coefficient_count(max_l))
    return;

  auto  index = longitude_index * output_resolution.x + latitude_index;
  auto& point = output_points[index];

  if (coefficient_index == 0)
  {
    point.y = 2 * M_PI * longitude_index / output_resolution.x;
    point.z =     M_PI * latitude_index  / output_resolution.y;
  }
  atomicAdd(&point.x, evaluate(coefficient_index, point.y, point.z) * coefficients[coefficient_index]);

  output_indices[index    ] =  longitude_index                            * output_resolution.y +  latitude_index,
  output_indices[index + 1] =  longitude_index                            * output_resolution.y + (latitude_index + 1) % output_resolution.y,
  output_indices[index + 2] = (longitude_index + 1) % output_resolution.x * output_resolution.y + (latitude_index + 1) % output_resolution.y,
  output_indices[index + 3] = (longitude_index + 1) % output_resolution.x * output_resolution.y +  latitude_index;
}
// Call on a dimensions.x x dimensions.y x dimensions.z 3D grid.
template<typename precision, typename point_type>
GLOBAL void sample_sums(
  const uint3        dimensions       ,
  const unsigned int max_l            ,
  const uint2        output_resolution,
  const precision*   coefficients     ,
  point_type*        output_points    ,
  unsigned int*      output_indices   )
{
  auto x = blockIdx.x * blockDim.x + threadIdx.x;
  auto y = blockIdx.y * blockDim.y + threadIdx.y;
  auto z = blockIdx.z * blockDim.z + threadIdx.z;
  
  if (x > dimensions.x || y > dimensions.y || z > dimensions.z)
    return;
  
  auto volume_index        = z + dimensions.z * (y + dimensions.y * x);
  auto coefficients_offset = volume_index * coefficient_count(max_l);
  auto samples_offset      = volume_index * output_resolution.x * output_resolution.y;

  sample_sum<<<dim3(output_resolution.x, output_resolution.y, coefficient_count(max_l)), 1>>>(
    max_l, 
    output_resolution, 
    coefficients   + coefficients_offset, 
    output_points  +     samples_offset ,
    output_indices + 4 * samples_offset );
}

// Call on a output_resolution.x x output_resolution.y 2D grid.
template<typename precision, typename point_type>
GLOBAL void sample(
  const unsigned int l                ,
  const int          m                ,
  const uint2        output_resolution,
  point_type*        output_points    ,
  unsigned int*      output_indices   )
{
  auto longitude_index = blockIdx.x * blockDim.x + threadIdx.x;
  auto latitude_index  = blockIdx.y * blockDim.y + threadIdx.y;
  
  if (longitude_index > output_resolution.x ||
      latitude_index  > output_resolution.y )
    return;

  auto  index = longitude_index * output_resolution.x + latitude_index;
  auto& point = output_points[index];
  
  point.y = 2 * M_PI * longitude_index / output_resolution.x;
  point.z =     M_PI * latitude_index  / output_resolution.y;
  point.x = evaluate(l, m, point.y, point.z);

  output_indices[index    ] =  longitude_index                            * output_resolution.y +  latitude_index,
  output_indices[index + 1] =  longitude_index                            * output_resolution.y + (latitude_index + 1) % output_resolution.y,
  output_indices[index + 2] = (longitude_index + 1) % output_resolution.x * output_resolution.y + (latitude_index + 1) % output_resolution.y,
  output_indices[index + 3] = (longitude_index + 1) % output_resolution.x * output_resolution.y +  latitude_index;
}
// Call on a dimensions.x x dimensions.y x dimensions.z 3D grid.
template<typename precision, typename point_type>
GLOBAL void sample(
  const uint3        dimensions       ,
  const unsigned int l                ,
  const int          m                ,
  const uint2        output_resolution,
  point_type*        output_points    ,
  unsigned int*      output_indices   )
{
  auto x = blockIdx.x * blockDim.x + threadIdx.x;
  auto y = blockIdx.y * blockDim.y + threadIdx.y;
  auto z = blockIdx.z * blockDim.z + threadIdx.z;
  
  if (x > dimensions.x || y > dimensions.y || z > dimensions.z)
    return;

  auto volume_index   = z + dimensions.z * (y + dimensions.y * x);
  auto samples_offset = volume_index * output_resolution.x * output_resolution.y;

  sample<<<dim3(output_resolution.x, output_resolution.y), 1>>>(
    l,
    m,
    output_resolution,
    output_points  +     samples_offset,
    output_indices + 4 * samples_offset);
}

// Call on a coefficient_count x coefficient_count x coefficient_count 3D grid.
// Based on Modern Quantum Mechanics 2nd Edition page 216 by Jun John Sakurai.
template<typename precision, typename atomics_precision = float>
GLOBAL void product(
  const unsigned int coefficient_count,
  const precision*   lhs_coefficients ,
  const precision*   rhs_coefficients ,
  atomics_precision* out_coefficients )
{
  auto lhs_index = blockIdx.x * blockDim.x + threadIdx.x;
  auto rhs_index = blockIdx.y * blockDim.y + threadIdx.y;
  auto out_index = blockIdx.z * blockDim.z + threadIdx.z;
  
  if (lhs_index > coefficient_count ||
      rhs_index > coefficient_count ||
      out_index > coefficient_count)
    return;

  auto lhs_lm   = coefficient_lm(lhs_index);
  auto rhs_lm   = coefficient_lm(rhs_index);
  auto out_lm   = coefficient_lm(out_index);
  auto cg1      = clebsch_gordan<atomics_precision>(lhs_lm.x, rhs_lm.x, out_lm.x, 0, 0, 0);
  auto cg2      = clebsch_gordan<atomics_precision>(lhs_lm.x, rhs_lm.x, out_lm.x, lhs_lm.y, rhs_lm.y, out_lm.y);
  auto coupling = sqrt((2 * lhs_lm.x + 1) * (2 * rhs_lm.x + 1) / (4 * M_PI * (2 * out_lm.x + 1))) * cg1 * cg2;

  atomicAdd(&out_coefficients[out_index], atomics_precision(coupling * lhs_coefficients[lhs_index] * rhs_coefficients[rhs_index]));
}
// Call on a dimensions.x x dimensions.y x dimensions.z 3D grid.
template<typename precision, typename atomics_precision = float>
GLOBAL void product(
  const uint3        dimensions       ,
  const unsigned int coefficient_count,
  const precision*   lhs_coefficients ,
  const precision*   rhs_coefficients ,
  atomics_precision* out_coefficients )
{
  auto x = blockIdx.x * blockDim.x + threadIdx.x;
  auto y = blockIdx.y * blockDim.y + threadIdx.y;
  auto z = blockIdx.z * blockDim.z + threadIdx.z;

  if (x > dimensions.x || y > dimensions.y || z > dimensions.z)
    return;

  auto volume_index        = z + dimensions.z * (y + dimensions.y * x);
  auto coefficients_offset = volume_index * coefficient_count;

  product<<<dim3(coefficient_count, coefficient_count, coefficient_count), 1>>>(
    coefficient_count,
    lhs_coefficients + coefficients_offset,
    rhs_coefficients + coefficients_offset,
    out_coefficients + coefficients_offset);
}
}

#endif
