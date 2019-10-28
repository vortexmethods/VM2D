link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

SUBDIRLIST(VMlibSubdir_VM2D ${CMAKE_SOURCE_DIR}/src/VM2D)
SUBDIRLIST(VMlibSubdir_VM3D ${CMAKE_SOURCE_DIR}/src/VM3D)
SUBDIRLIST(VMlibSubdir_VMlib ${CMAKE_CURRENT_SOURCE_DIR})
LIST(APPEND VMlibDirs ${VMlibSubdir_VMlib} ${VMlibSubdir_VM2D} ${VMlibSubdir_VM3D})

#message(STATUS "VMlib_Dirs: ${VMlibDirs}")

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${VMlibDirs})

FILE(GLOB VMlib_ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/*/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cuh
  )

#MPI
if(WIN32)
   set(MPI_GUESS_LIBRARY_NAME "MSMPI" CACHE STRING "Override MPI library name, e.g., MSMPI, MPICH2")
endif()
   find_package(MPI REQUIRED)   
   include_directories(${MPI_CXX_INCLUDE_PATH})
   link_libraries(${MPI_CXX_LIBRARIES})


#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

IF(EXISTS "${CMAKE_SOURCE_DIR}/src/VM2D/CMakeLists.txt")
	#message(STATUS "-------------------------2D CODE-------------------------")
	add_library(VMlib2 STATIC ${VMlib_ALLFILES})
    target_compile_definitions(VMlib2 PUBLIC -DCODE2D)	
ENDIF()

IF(EXISTS "${CMAKE_SOURCE_DIR}/src/VM3D/CMakeLists.txt")
	#message(STATUS "-------------------------3D CODE-------------------------")
	add_library(VMlib3 STATIC ${VMlib_ALLFILES})
    target_compile_definitions(VMlib3 PUBLIC -DCODE3D)	
ENDIF()  
  

