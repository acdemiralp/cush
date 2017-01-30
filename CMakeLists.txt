cmake_minimum_required(VERSION 3.0.2)

# Setup project metadata.
set(ProjectName "cush")
project(${ProjectName} VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CUDA_NVCC_FLAGS "-gencode=arch=compute_50,code=sm_50 ${CUDA_NVCC_FLAGS}")

# Setup dependencies.
include("cmake/_dependencies.cmake")

# Setup sources.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include("cmake/_sources.cmake")

# Create source groups.
include("cmake/assign_source_group.cmake")
assign_source_group(${ProjectSources})

# Specify library.
cuda_add_library(${ProjectName} ${ProjectSources})

# Link libraries.
target_link_libraries(${ProjectName} ${ProjectLibraries})

# Add tests.
add_subdirectory(./tests)