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
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
message("${sdsl-lite_INCLUDE}")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# bsort (header-only, so here we just check if we can build it)
ExternalProject_Add(bsort
  GIT_REPOSITORY "https://github.com/ekg/bsort.git"
  GIT_TAG "e80bdd183ecca23e75901d2bdabe9920b78f1590"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(bsort SOURCE_DIR)
set(bsort_INCLUDE "${SOURCE_DIR}/src")
set(bsort_LIB "${SOURCE_DIR}/lib")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# gzip file reader, header only
ExternalProject_Add(gzipreader
  GIT_REPOSITORY "https://github.com/gatoravi/gzip_reader.git"
  GIT_TAG "0ef26c0399e926087f9d6c4a56067a7bf1fc4f5e"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(gzipreader SOURCE_DIR)
set(gzipreader_INCLUDE "${SOURCE_DIR}")

set(CMAKE_BUILD_TYPE Debug) 

# set up our target executable and specify its dependencies and includes
add_executable(seqwish
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  ${CMAKE_SOURCE_DIR}/src/seqindex.cpp
  ${CMAKE_SOURCE_DIR}/src/paf.cpp
  ${CMAKE_SOURCE_DIR}/src/cigar.cpp
  ${CMAKE_SOURCE_DIR}/src/alignments.cpp
  ${CMAKE_SOURCE_DIR}/src/pos.cpp
  ${CMAKE_SOURCE_DIR}/src/transclosure.cpp
  ${CMAKE_SOURCE_DIR}/src/links.cpp
  ${CMAKE_SOURCE_DIR}/src/compact.cpp
  ${CMAKE_SOURCE_DIR}/src/dna.cpp
  ${CMAKE_SOURCE_DIR}/src/gfa.cpp
  )
add_dependencies(seqwish bsort)
add_dependencies(seqwish tayweeargs)
add_dependencies(seqwish sdsl-lite)
add_dependencies(seqwish gzipreader)
target_include_directories(seqwish PUBLIC
  "${bsort_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${gzipreader_INCLUDE}")
target_link_libraries(seqwish
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${bsort_LIB}/libbsort.a"
  z)
set(CMAKE_EXE_LINKER_FLAGS "-static")
