###################################################################################
#
# Vessel tree construction
#
# Copyright (c) 2020 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

cmake_minimum_required(VERSION 3.16)

option(enable_OPENMP "Enable OpenMP" ON)
option(enable_GLOG "Enable GLOG" OFF)

#-----compiler options
if(NOT DEFINED compiler)
  set(base_compile_options "-Wall -Wextra -O3")
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -fopenmp")
   endif()
endif()

if(compiler STREQUAL intel)
  set(CMAKE_C_COMPILER "icc")
  set(CMAKE_CXX_COMPILER "icpc")
  set(base_compile_options "-Wall -Wextra -fast")
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -qopenmp")
  endif()
elseif(DEFINED compiler)
  message(FATAL ERROR "invalid compiler")
endif()
string(REPLACE " " ";" compile_options ${base_compile_options})
#-----------------------

project(WSS_CSF C CXX)
set(CMAKE_CXX_STANDARD 11)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

#change default install prefix
set(INSTALL_DIR "/usr/local/CoRotationalShell")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH "..." FORCE)
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

#######
# Display options
#######
message(" ")
message("${PROJECT_NAME}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "COMPILE_OPTIONS       : "      ${compile_options})
message(" ")

# HDF5
find_package (HDF5 REQUIRED COMPONENTS C CXX)
IF(HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  MESSAGE("find HDF5")
  MESSAGE(${HDF5_INCLUDE_DIRS})
  MESSAGE(${HDF5_LIBRARIES})
  MESSAGE(${HDF5_CXX_LIBRARIES})
ELSE()
  MESSAGE(FATAL_ERROR "not find HDF5")
ENDIF(HDF5_FOUND)

#######
# Sub Directories
#######

# add_subdirectory(include)
add_subdirectory(src)

set(main_target "Main")

add_executable(${main_target} main.cpp)
target_sources(${main_target}
  PRIVATE
    main.cpp
)

target_include_directories(${main_target}
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${main_target} PRIVATE
  ${HDF5_LIBRARIES}
  ${HDF5_CXX_LIBRARIES}
  Base
  H5
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})

