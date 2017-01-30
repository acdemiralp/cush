# Register bundled modules.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include Cuda and CuSolver.
find_package(Cuda REQUIRED)
set(ProjectLibraries ${ProjectLibraries} "${CUDA_CUBLAS_LIBRARIES};${CUDA_cusolver_LIBRARY};${CUDA_cudadevrt_LIBRARY}")