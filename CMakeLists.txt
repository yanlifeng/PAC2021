# @file CMakeLists.txt
# @author Fande Yu 
# @version 0.1 
# @copyright 
#
# ChangeLog
# AUTHOR      | TIME       | VERSION       | DESCRIPTION
# ------      | ----       | -------       | -----------
# Fande Yu    | 2020/10/14 | 0.1.20201014  | Create this file 
# Fande Yu	  | 2020/12/26 | 0.1.20201226  | Build Seperatedly

#Notice
#1. export CC=gcc and CXX=g++ ahead of running this CMake file.
#2. If this program will be run on multi-nodes where NCCL may be necessary, this file need to be modified.

# CMAKE REQUIREMENT
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

#set(CMAKE_C_COMPILER "/usr/local/Cellar/gcc@8/8.4.0_1/bin/gcc-8")
#set(CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc@8/8.4.0_1/bin/g++-8")

option(USE_INTEL_COMPILER "Whether to use intel compiler?" ON)
option(USE_INTEL_FFTW "Whether to use intel fftw?" ON)

if (USE_INTEL_COMPILER)
    set(CMAKE_C_COMPILER icx)
    set(CMAKE_CXX_COMPILER icpx)
endif ()

project(TomoProject)

set(TomoProject_VERSION_MAJOR 0)
set(TomoProject_VERSION_MINOR 1)
set(TomoProject_VERSION_ADDIT 20201014)

#Get git commit
execute_process(
        COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE COMMIT_VERSION
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Welcome to TomoProject ${TomoProject_VERSION_MAJOR}.${TomoProject_VERSION_MINOR}.${TomoProject_VERSION_ADDIT}!")


#For TomoConfig.h.in
set(COMMIT_VERSION_QUOTE "\"${COMMIT_VERSION}\"")


#set(CMAKE_C_COMPILER gcc)
#set(CMAKE_CXX_COMPILER g++)
set(COMMON_FLAGS "-std=c++11 ${COMMON_FLAGS}")

#For FFTW
option(SINGLE_PRECISION "Whether to use single float precision in calculation?" ON)

#Only option is useful for #cmakedefine in TomoConfig.h.in
option(CPU_VERSION "Whether to use CPU to calculate?" ON)


#Add compilation flag part1: build type
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "release")
endif ("${CMAKE_BUILD_TYPE}" STREQUAL "")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "release")
    set(COMMON_FLAGS "${COMMON_FLAGS} -fopenmp -O2")
    message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}, compile ${CMAKE_PROJECT_NAME} with ${COMMON_FLAGS} flag.")
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "debug")
    #TODO 拆成additional_flags
    set(COMMON_FLAGS "${COMMON_FLAGS} -fopenmp -O2")
    set(COMMON_FLAGS "${COMMON_FLAGS} -g")
    message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}, compile ${CMAKE_PROJECT_NAME} with ${COMMON_FLAGS} flag.")
endif ("${CMAKE_BUILD_TYPE}" STREQUAL "release")

if ("${CPU_VERSION}")
    message(STATUS "Try to build Tomo-Project in CPU version.")
endif ("${CPU_VERSION}")

if (USE_INTEL_COMPILER)
    set(COMMON_FLAGS "${COMMON_FLAGS} -fp-model precise")
endif ()
set(COMMON_FLAGS "${COMMON_FLAGS} -lz -lpthread -lm ")

set(CMAKE_C_FLAGS "${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS "${COMMON_FLAGS}")

message(STATUS "CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")

if (USE_INTEL_FFTW)
    set(FFTW_PATH "/home/PAC20217111/intel/oneapi/mkl/2021.3.0")
else ()
    set(FFTW_PATH "/home/PAC20217111/ylf/fftwTest/fftw")
endif ()
message(STATUS "FFTW Path: ${FFTW_PATH}")

include_directories(${FFTW_PATH}/include)

if (USE_INTEL_FFTW)
    include_directories(${FFTW_PATH}/include/fftw)
    link_directories(${FFTW_PATH}/lib/intel64)
else ()
    link_directories(${FFTW_PATH}/lib)
endif ()

include_directories(${CMAKE_SOURCE_DIR}/external/Eigen)
include_directories(${CMAKE_SOURCE_DIR}/src/include)


if ("${CPU_VERSION}")


    file(GLOB_RECURSE CORE_SOURCES ${CMAKE_SOURCE_DIR}/src/core/*.cpp)
    add_library(core_CPU STATIC ${CORE_SOURCES})
    if (USE_INTEL_COMPILER)
        if (USE_INTEL_FFTW)
            target_link_libraries(core_CPU mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 dl)
        else ()
            target_link_libraries(core_CPU fftw3f fftw3f_threads)
        endif ()
    else ()
        if (USE_INTEL_FFTW)
            target_link_libraries(core_CPU mkl_intel_lp64 mkl_gnu_thread mkl_core gomp dl)
        else ()
            target_link_libraries(core_CPU fftw3f fftw3f_threads)
        endif ()
    endif ()
    file(GLOB_RECURSE UTIL_SOURCES ${CMAKE_SOURCE_DIR}/src/util/*.cpp)
    add_library(util STATIC ${UTIL_SOURCES})

    file(GLOB_RECURSE MAIN_SOURCES ${CMAKE_SOURCE_DIR}/src/main/*.cpp)

    foreach (MAINSRC ${MAIN_SOURCES})
        get_filename_component(MAINNAME ${MAINSRC} NAME_WE) # file name without directories or longest extension
        set(MAINNAME "${MAINNAME}_cpu")
        add_executable(${MAINNAME} ${MAINSRC})
        if (USE_INTEL_COMPILER)
            if (USE_INTEL_FFTW)
                target_link_libraries(${MAINNAME} core_CPU util mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 dl)
            else ()
                target_link_libraries(${MAINNAME} core_CPU util fftw3f fftw3f_threads)
            endif ()
        else ()
            if (USE_INTEL_FFTW)
                target_link_libraries(${MAINNAME} core_CPU util mkl_intel_lp64 mkl_gnu_thread mkl_core gomp dl)
            else ()
                target_link_libraries(${MAINNAME} core_CPU util fftw3f fftw3f_threads)
            endif ()
        endif ()
    endforeach ()

endif ("${CPU_VERSION}")
