cmake_minimum_required(VERSION 3.12)

# Guard: when included via add_subdirectory() from the parent, skip the
# standalone project() call so we don't clobber its settings.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    #project(FMM VERSION 1.0.0 LANGUAGES CXX CUDA)
    project(FMM VERSION 1.0.0 LANGUAGES CXX)
    # Only force Release build when built standalone
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
endif()

find_package(OpenMP REQUIRED)
find_package(TBB REQUIRED)

# MPI is optional — enable by passing -DFMM_USE_MPI=ON
option(FMM_USE_MPI "Enable MPI support in fmmlib" OFF)
if(FMM_USE_MPI)
    find_package(MPI REQUIRED)
endif()

FILE(GLOB fmmlib_src
    ${CMAKE_CURRENT_SOURCE_DIR}/cpu/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cpu/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/common/*.h)

add_library(fmmlib ${fmmlib_src})

target_compile_features(fmmlib PRIVATE cxx_std_20)
target_compile_options(fmmlib PUBLIC -mavx2)
# Suppress warnings from this external contribution
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fmmlib PRIVATE
        -w  # Suppress all warnings
    )
endif()
target_include_directories(fmmlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(FMM_USE_MPI)
    target_compile_definitions(fmmlib PUBLIC FMM_MPI)
    target_link_libraries(fmmlib PUBLIC TBB::tbb TBB::tbbmalloc OpenMP::OpenMP_CXX MPI::MPI_CXX)
else()
    target_link_libraries(fmmlib PUBLIC TBB::tbb TBB::tbbmalloc OpenMP::OpenMP_CXX)
endif()

# Standalone example executable — only built when this is the top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_executable(FMM ${CMAKE_CURRENT_SOURCE_DIR}/src/examples/main.cpp)
    target_link_libraries(FMM PUBLIC fmmlib)
endif()
