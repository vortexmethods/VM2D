link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

SUBDIRLIST(VMlibSubdir_VM2D ${CMAKE_SOURCE_DIR}/src/VM2D)
SUBDIRLIST(VMlibSubdir_VM3D ${CMAKE_SOURCE_DIR}/src/VM3D)
SUBDIRLIST(VMlibSubdir_VMlib ${CMAKE_CURRENT_SOURCE_DIR})
LIST(APPEND VMlibDirs ${VMlibSubdir_VMlib} ${VMlibSubdir_VM2D} ${VMlibSubdir_VM3D})

message(STATUS "VMlib_Dirs: ${VMlibDirs}")

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${VMlibDirs})

FILE(GLOB VMlib_ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/*/*.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cuh
  )

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
   link_libraries(${MPI_CXX_LIBRARIES})
else()
   find_package(MPI REQUIRED)
   include_directories(${MPI_CXX_INCLUDE_PATH})
   link_libraries(${MPI_CXX_LIBRARIES})
endif()


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
  

