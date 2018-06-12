# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# We build using c++11
set(CMAKE_CXX_STANDARD 11)

# Project's name
project(seqwish)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "d52aa9a71513d132e30c09491b5899af449ebb94"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
message("${sdsl-lite_INCLUDE}")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")

# bsort (header-only, so here we just check if we can build it)
ExternalProject_Add(bsort
  GIT_REPOSITORY "https://github.com/ekg/bsort.git"
  GIT_TAG "0e4076e176d8777d6287bdbe373078415f858903"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(bsort SOURCE_DIR)
set(bsort_INCLUDE "${SOURCE_DIR}/src")

# TODO link into bsort and sdsl-lite to verify build works
add_executable(seqwish ${CMAKE_SOURCE_DIR}/src/main.cpp)
add_dependencies(seqwish bsort)
add_dependencies(seqwish sdsl-lite)
target_include_directories(seqwish PUBLIC
  "${bsort_INCLUDE}"
  "${sdsl-lite_INCLUDE}")
target_link_libraries(seqwish "${sdsl-lite_LIB}/libsdsl.a")
set(CMAKE_EXE_LINKER_FLAGS "-static")
