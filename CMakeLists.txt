# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# We build using c++11
set(CMAKE_CXX_STANDARD 11)

# Project's name
project(seqwish)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
#ExternalProject_Add(sdsl-lite
#  SOURCE_DIR ${CMAKE_SOURCE_DIR}/deps/sdsl-lite)
#set(SDSL_INSTALL_PREFIX
#CMAKE_INSTALL_PREFIX
ExternalProject_Add(sdsl-lite
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "d52aa9a71513d132e30c09491b5899af449ebb94"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")

ExternalProject_Add(bsort
  GIT_REPOSITORY "https://github.com/ekg/bsort.git"
  GIT_TAG "0e4076e176d8777d6287bdbe373078415f858903"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")

add_executable(seqwish ${CMAKE_SOURCE_DIR}/src/main.cpp)
