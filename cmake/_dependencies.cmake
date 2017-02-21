# Register bundled modules.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include Cuda.
find_package(Cuda REQUIRED)
set(ProjectLibraries ${ProjectLibraries} "${CUDA_cudadevrt_LIBRARY}")