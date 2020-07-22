link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

SUBDIRLIST(VM2Dsubdir_VM2D ${CMAKE_CURRENT_SOURCE_DIR})
SUBDIRLIST(VM2Dsubdir_VMlib ${CMAKE_SOURCE_DIR}/src/VMlib)
SUBDIRLIST(VM2Dsubdir_VMcuda ${CMAKE_SOURCE_DIR}/src/VMcuda)
LIST(APPEND VM2Ddirs ${VM2Dsubdir_VM2D} ${VM2Dsubdir_VMlib} ${VM2Dsubdir_VMcuda})

#message(STATUS "VM2D_dirs: ${VM2Ddirs}")

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${VM2Ddirs})

FILE(GLOB VM2D_ALLFILES
      ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/*/*.h
    )

#MPI
if(WIN32)
   if (${CMAKE_VS_PLATFORM_TOOLSET} MATCHES "Intel")
      #message(STATUS "INTEL COMPILER")	  
   else()
      set(MPI_GUESS_LIBRARY_NAME "MSMPI" CACHE STRING "Override MPI library name, e.g., MSMPI, MPICH2")
	  #message(STATUS "NOT INTEL COMPILER")	  
   endif()  
endif()
   find_package(MPI REQUIRED)   
   include_directories(${MPI_CXX_INCLUDE_PATH})
   link_libraries(${MPI_CXX_LIBRARIES})
 
#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	
	#message(STATUS "OMP directory ${OpenMP_LIBRARY_DIR}")

endif()

 
add_executable(VM2D ${VM2D_ALLFILES})

target_link_libraries(VM2D debug VMlib2)
target_link_libraries(VM2D optimized VMlib2)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(WIN32)
        target_link_libraries(VM2D debug "c:/Program Files/LLVM/lib/libomp.lib")
		target_link_libraries(VM2D optimized "c:/Program Files/LLVM/lib/libomp.lib")
	endif()
endif()	

#CUDA
if(CUDA_FOUND OR CMAKE_CUDA_COMPILER)
   target_link_libraries(VM2D debug VMcuda)
   target_link_libraries(VM2D optimized VMcuda)
endif()
