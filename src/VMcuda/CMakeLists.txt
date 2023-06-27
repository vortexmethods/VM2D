link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

SUBDIRLIST(VMcudasubdir_VMcuda ${CMAKE_CURRENT_SOURCE_DIR})
SUBDIRLIST(VMcudasubdir_VMlib ${CMAKE_SOURCE_DIR}/src/VMlib)
LIST(APPEND VMcudadirs ${VMcudasubdir_VMcuda} ${VMcudasubdir_VMlib})

#message(STATUS "VMcuda_dirs: ${VMcudadirs}")

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${VMcudadirs})

include_directories(${CUDA_INCLUDE_DIRS})

FILE(GLOB VMcuda_ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cuh
    )

#CUDA
if (CMAKE_CUDA_COMPILER)
    add_library(VMcuda ${VMcuda_ALLFILES})
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
endif()