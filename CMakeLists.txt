cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(TIOGA CXX)

option(BUILD_SHARED_LIBS "Build shared libraries (default: off)" ON)
option(BUILD_TIOGA_EXE "Build tioga driver code (default: off)" OFF)
option(BUILD_GRIDGEN_EXE "Build grid generator code (default: off)" OFF)
option(TIOGA_HAS_NODEGID "Support node global IDs (default: on)" ON)
option(TIOGA_ENABLE_TIMERS "Track timing information for TIOGA (default: off)" OFF)
option(TIOGA_OUTPUT_STATS "Output statistics for TIOGA holecutting (default: off)" OFF)
option(TIOGA_ENABLE_CUDA "Enable CUDA support (default: off)" OFF)
option(TIOGA_ENABLE_HIP "Enable AMD HIP support (default: off)" OFF)
option(TIOGA_ENABLE_ARBORX "Enable ArborX (default: off)" OFF)

# CUDA specific options
#set(TIOGA_CUDA_SM "70" CACHE STRING "CUDA arch option")
#set(TIOGA_HIP_ARCH "gfx906,gfx908" CACHE STRING "HIP arch option")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)

if (TIOGA_ENABLE_ARBORX)
  add_definitions(-DTIOGA_USE_ARBORX)
  find_package(ArborX REQUIRED)
endif()

#if (TIOGA_ENABLE_CUDA)
#  enable_language(CUDA)
#  if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "9.0")
#    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION}. "
#      "Please use CUDA toolkit version 9.0 or newer")
#  endif()
#  set(TIOGA_CUDA_FLAGS "--Wno-deprecated-gpu-targets -m64")
#  set(TIOGA_CUDA_ARCH_FLAGS "-gencode=arch=compute_${TIOGA_CUDA_SM},code=sm_${TIOGA_CUDA_SM} -expt-extended-lambda")
#  set(CMAKE_CUDA_STANDARD 11)
#  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
#  set(CMAKE_CUDA_EXTENSIONS OFF)
#endif()

#if (TIOGA_ENABLE_HIP)
#  if (NOT DEFINED HIP_PATH)
#    set(HIP_PATH $ENV{HIP_PATH})
#  endif()
#  set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" "${CMAKE_MODULE_PATH}")
#  find_package(HIP REQUIRED)
#  if ("${HIP_COMPILER}" STREQUAL "hcc")
#    message(FATAL_ERROR "TIOGA does not support deprecated HCC compiler. Use newer ROCm")
#  endif()
#  message(STATUS "HIP: ${HIP_VERSION}, platform = ${HIP_PLATFORM}, compiler = ${HIP_COMPILER}")
#  set(TIOGA_HIP_ARCH_FLAGS "-m64 --amdgpu-target=${TIOGA_HIP_ARCH}")
#endif()

add_subdirectory(src)

# Optionally build driver exe and gridGen if the user requests it
if (BUILD_TIOGA_EXE)
  add_subdirectory(driver)
endif()

if (BUILD_GRIDGEN_EXE)
  add_subdirectory(gridGen)
endif()

# CMake installation configuration

install(EXPORT TIOGALibraries
  DESTINATION lib/cmake/TIOGA
  FILE TIOGALibraries.cmake)

# Create TIOGA config so that other codes can find TIOGA
include(CMakePackageConfigHelpers)
set(INCLUDE_INSTALL_DIR include/)
set(LIB_INSTALL_DIR lib/)
configure_package_config_file(
  cmake/TIOGAConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/TIOGAConfig.cmake
  INSTALL_DESTINATION lib/cmake/TIOGA
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/TIOGAConfig.cmake
  DESTINATION lib/cmake/TIOGA)
