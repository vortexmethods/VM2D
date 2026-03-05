include(FindPackageHandleStandardArgs)

find_path(cuDNN_INCLUDE_DIR
    NAMES cudnn.h
    PATHS ${CMAKE_INSTALL_PREFIX} ${CMAKE_SOURCE_DIR} ${CUDAToolkit_INCLUDE_DIRS} ${CUDAToolkit_ROOT}
    PATH_SUFFIXES include cuda/include
    REQUIRED
)
find_path(cuDNN_LIBRARY_DIR
    NAMES
        "libcudnn.so"
        "libcudnn.so.${MAJOR_NUMBER}"
        "libcudnn_graph_static.a"
        "cudnn.lib"
        "cudnn64_${MAJOR_NUMBER}.dll"
        "cudnn_graph_static64_${MAJOR_NUMBER}.lib"
    PATHS ${CMAKE_INSTALL_PREFIX} ${CUDAToolkit_LIBRARY_DIR} ${CUDAToolkit_ROOT}
    PATH_SUFFIXES lib lib64 lib/x64
    REQUIRED
)



message(STATUS "Found cuDNN: ${cuDNN_LIBRARY_DIR}")

set(_shared_libs
    cudnn
    cudnn_ops
    cudnn_adv
    cudnn_cnn
    cudnn_graph
    cudnn_heuristic
    cudnn_engines_runtime_compiled
    cudnn_engines_precompiled
)
set(_static_libs
    cudnn_ops_static
    cudnn_adv_static
    cudnn_cnn_static
    cudnn_graph_static
    cudnn_heuristic_static
    cudnn_engines_runtime_compiled_static
    cudnn_engines_precompiled_static
)

foreach(shared IN LISTS _shared_libs)
    find_library(${shared}_loc
        NAMES ${shared}
        PATHS ${cuDNN_LIBRARY_DIR}
        NO_DEFAULT_PATH
    )

    add_library(${shared} SHARED IMPORTED)
    target_include_directories(${shared} INTERFACE ${cuDNN_INCLUDE_DIR})
    if(NOT MSVC)
        set_property(TARGET ${shared} PROPERTY
            IMPORTED_LOCATION "${${shared}_loc}"
        )
    else()
        set_property(TARGET ${shared} PROPERTY
            IMPORTED_IMPLIB "${${shared}_loc}"
        )
    endif()
endforeach()

if(NOT WIN32)
    set_property(TARGET cudnn_graph PROPERTY INTERFACE_LINK_LIBRARIES ZLIB::ZLIB)
endif()

if(NOT MSVC OR NOT BUILD_WIN_SHARED_LIBS)
    # find paths to static libs and create imports
    foreach(static IN LISTS _static_libs)
        find_library(${static}_loc
            NAMES ${static} ${static}64_${MAJOR_NUMBER}
            PATHS ${cuDNN_LIBRARY_DIR}
            NO_DEFAULT_PATH
        )
        if(${static}_loc)
            add_library(${static} STATIC IMPORTED)
            target_include_directories(${static} INTERFACE ${cuDNN_INCLUDE_DIR})
            set_property(TARGET ${static} PROPERTY
                        IMPORTED_LOCATION "${${static}_loc}"
            )
        endif()
    endforeach()

    #set_property(TARGET cudnn_graph_static
    #    PROPERTY INTERFACE_LINK_LIBRARIES
    #    ZLIB::ZLIB
    #)
	
    #if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 11.5)
    #    set_property(TARGET cudnn_graph_static
    #        APPEND PROPERTY INTERFACE_LINK_LIBRARIES
    #        CUDA::nvrtc_static
    #        CUDA::nvrtc_builtins_static
    #        CUDA::nvptxcompiler_static
    #    )
    #endif()
endif()
