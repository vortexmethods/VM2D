#cmake .. -DPATH_TO_EIGEN="C:/Program files (x86)/Eigen3/include/eigen3"

cmake_minimum_required(VERSION 2.8)
project(VM2D)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

set(SRC_LIST src/VM2D/VM2D.cpp)

FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h
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

#MESSAGE( STATUS "variable is: " ${dirlist} )

include_directories(${dirlist})

#EIGEN
include_directories(${PATH_TO_EIGEN})

#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()	

