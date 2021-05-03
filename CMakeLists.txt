cmake_minimum_required(VERSION 3.13.0)
project(VM)

set(CudaOldIncludeMethod FALSE)

set(VM2D_DIR "build")


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")

set(VM_WORK_DIR "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")

message(STATUS "Work directory ${VM_WORK_DIR}")

#message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}")
#message(STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}")
#message(STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}")


include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11")
endif()

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


#CUDA OLD
if(CudaOldIncludeMethod)
   find_package(CUDA)   
#CUDA NEW
else()
   include(CheckLanguage)
   check_language(CUDA)
   if (CUDA_FOUND OR CMAKE_CUDA_COMPILER)
		enable_language(CUDA)   
   endif()
endif()

#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_53,code=sm_53 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_62,code=sm_62")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_53,code=sm_53 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_62,code=sm_62 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_72,code=sm_72 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80")


#COMMENT THIS LINE IF THERE IS INSTALLED CUDA SDK BUT YOU DON'T WANT TO USE GPU
if (CUDA_FOUND OR CMAKE_CUDA_COMPILER)
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
add_subdirectory(src/VMcuda)

IF(EXISTS "${CMAKE_SOURCE_DIR}/src/VM2D/CMakeLists.txt")
	add_subdirectory(src/VM2D)
	set_target_properties(VM2D PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
	set_target_properties(VM2D PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	set_target_properties(VM2D PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/${VM2D_DIR}")
		
	if (CUDA_FOUND OR CMAKE_CUDA_COMPILER)
		add_dependencies(VM2D VMlib2 VMcuda)
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
