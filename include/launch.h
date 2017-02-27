#ifndef CUSH_LAUNCH_H_
#define CUSH_LAUNCH_H_

#include <math.h>

#include <device_launch_parameters.h>

#include <decorators.h>

namespace cush
{
INLINE COMMON unsigned block_size_1d()
{
  return 32;
}
INLINE COMMON dim3     block_size_2d()
{
  return {16, 16, 1};
}
INLINE COMMON dim3     block_size_3d()
{
  return {8, 8, 8};
}

INLINE COMMON unsigned grid_size_1d(unsigned target_dimension )
{
  return ceil(float(target_dimension) / block_size_1d());
}
INLINE COMMON dim3     grid_size_2d(dim3     target_dimensions)
{
  auto block_size = block_size_2d();
  return {
    ceil(float(target_dimensions.x) / block_size.x),
    ceil(float(target_dimensions.y) / block_size.y),
    1 
  };
}
INLINE COMMON dim3     grid_size_3d(dim3     target_dimensions)
{
  auto block_size = block_size_3d();
  return{
    ceil(float(target_dimensions.x) / block_size.x),
    ceil(float(target_dimensions.y) / block_size.y),
    ceil(float(target_dimensions.z) / block_size.z)
  };
}
}

#endif
