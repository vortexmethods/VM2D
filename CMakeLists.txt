cmake_minimum_required(VERSION 2.8)
project(VM2D)

enable_language(CUDA)
set (CMAKE_CUDA_FLAGS "-arch=sm_35")
#set (CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

set(CMAKE_CONFIGURATION_TYPE "Debug;Release")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")	
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

set(SRC_LIST src/VM2D/VM2D.cpp)

FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cuh
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cuh
  )

add_executable(${PROJECT_NAME} ${ALLFILES})

#FOREACH(subdir ${CMAKE_CURRENT_SOURCE_DIR}/src)
#%  ADD_SUBDIRECTORY(${subdir})
#ENDFOREACH()

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()


SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)

#MESSAGE( STATUS "variable is: " ${CMAKE_CONFIGURATION_TYPE} )

include_directories(${dirlist})

include_directories(/opt/openmpi/intel/eth/include)

#EIGEN
#if(WIN32)
#set(PATH_TO_EIGEN "C:/Program files (x86)/Eigen3/include/eigen3")
#include_directories(${PATH_TO_EIGEN})
#else()
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
#find_package(Eigen3 REQUIRED)
set(EIGEN3_INCLUDE_DIR "./include")
include_directories(${EIGEN3_INCLUDE_DIR})
#endif()

#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()




# Require MPI for this project:
if(WIN32)
   #This is for finding MS-MPI.
   #set(_MPI_PREFIX_PATH)
   #list(APPEND _MPI_PREFIX_PATH "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\MPI;Path]/..")
   set(MPIEXEC "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe")
   # For building MPI programs the selected Visual Studio compiler is used, namely cl.exe.
   # So there is no need to set a specific MPI compiler.
   #set(MPI_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
   set(MPI_CXX_INCLUDE_PATH "C:/Program Files (x86)/Microsoft SDKs/MPI/Include")
   # Make sure the correct libraries (64-bit or 32-bit) are selected.
   # Decide between 32-bit and 64-bit libraries for Microsoft's MPI
   if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
     set(MS_MPI_ARCH_DIR x64)
   else()
     set(MS_MPI_ARCH_DIR x86)
   endif()
   set(MPI_CXX_LIBRARIES "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/${MS_MPI_ARCH_DIR}/msmpi.lib")
   set(MPI_C_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
   set(MPI_C_LIBRARIES "{${MPI_CXX_LIBRARIES}")
   set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
  
   include_directories(${MPI_CXX_INCLUDE_PATH})
   target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
else()
    find_package(MPI REQUIRED)
endif()

#MPI
#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})
#target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()	

