cmake_minimum_required(VERSION 3.13.0)
project(VM)

set(VM2D_DIR "build")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")

set(VM_WORK_DIR "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")

message(STATUS "Work directory ${VM_WORK_DIR}")

#message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}")
#message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}")
#message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}")


#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#if(COMPILER_SUPPORTS_CXX11)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11")
#endif()

message (STATUS "C++ compiler is: " ${CMAKE_CXX_COMPILER})

if (CMAKE_COMPILER_IS_INTEL OR (${CMAKE_VS_PLATFORM_TOOLSET} MATCHES "Intel") OR (${CMAKE_CXX_COMPILER} MATCHES "Intel")  OR (${CMAKE_CXX_COMPILER} MATCHES "intel") )
    message (STATUS "Flags for Intel C++ compiler")
	if (WIN32)
	   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /fp:except /fp:precise")
	else() 
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fp-model=except -fp-model=precise")
    endif()	   
else()
    message (STATUS "Flags for non-Intel C++ compiler")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()


#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


#MPI
if (WIN32)
   if (${CMAKE_VS_PLATFORM_TOOLSET} MATCHES "Intel")
      #message(STATUS "INTEL COMPILER")
	  #message(STATUS "CXX: " ${CXX})
   else()
      set(MPI_GUESS_LIBRARY_NAME "MSMPI" CACHE STRING "Override MPI library name, e.g., MSMPI, MPICH2")
	  #message(STATUS "NOT INTEL COMPILER")
	  #message(STATUS "CXX: " ${CXX})
   endif()
endif()
  
find_package(MPI)
if (MPI_FOUND)
	include_directories(${MPI_CXX_INCLUDE_PATH})
	link_libraries(${MPI_CXX_LIBRARIES})
	add_definitions(-DUSE_MPI)
endif()


include(CheckLanguage)
check_language(CUDA)
if (CUDA_FOUND OR CMAKE_CUDA_COMPILER)
   enable_language(CUDA)   
   #set (ARCHITECTURE "60")
   set (ARCHITECTURE "60;70")
   #set (ARCHITECTURE "70;80")  
   
   #COMMENT THIS LINE IF THERE IS INSTALLED CUDA SDK BUT YOU DON'T WANT TO USE GPU   
   add_definitions(-DUSE_CUDA)   
endif()

if(MSVC)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

IF(EXISTS "${CMAKE_SOURCE_DIR}/src/VM2D/CMakeLists.txt")
	message(STATUS "-------------------------2D CODE-------------------------")
ENDIF()

IF(EXISTS "${CMAKE_SOURCE_DIR}/src/VM3D/CMakeLists.txt")
	message(STATUS "-------------------------3D CODE-------------------------")
ENDIF()


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

add_subdirectory(src/VMlib)
if (CUDA_FOUND OR CMAKE_CUDA_COMPILER)
    add_subdirectory(src/VMcuda)
endif() 

IF(EXISTS "${CMAKE_SOURCE_DIR}/src/VM2D/CMakeLists.txt")
	add_subdirectory(src/VM2D)
	set_target_properties(VM2D PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
	set_target_properties(VM2D PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	set_target_properties(VM2D PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")
		
	if (CUDA_FOUND OR CMAKE_CUDA_COMPILER)		
		set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
		set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
		
		message(STATUS "-------------------------CUDA FLAGS-------------------------")
		
		add_dependencies(VM2D VMlib2 VMcuda)
		set_target_properties(VMcuda PROPERTIES CUDA_ARCHITECTURES "${ARCHITECTURE}")
		set_target_properties(VM2D PROPERTIES CUDA_ARCHITECTURES "${ARCHITECTURE}")
		target_link_libraries(VM2D cusolver)
	else()
		add_dependencies(VM2D VMlib2)
	endif()
	
	set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT VM2D)
ENDIF()

IF(EXISTS "${CMAKE_SOURCE_DIR}/src/VM3D/CMakeLists.txt")
	add_subdirectory(src/VM3D)
	set_target_properties(VM3D PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
	set_target_properties(VM3D PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	set_target_properties(VM3D PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")

	add_dependencies(VM3D VMlib3)
	#set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT VM3D)
ENDIF()
