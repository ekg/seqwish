# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(seqwish)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -mcx16 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -mcx16 -g")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # assumes clang build
  # we can't reliably detect when we're using clang, so for the time being we assume
  # TODO: can't we though?
  # adapted from https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
  # find_package(OpenMP) does not work reliably on macOS, so we do its work ourselves
  set (OpenMP_C "${CMAKE_C_COMPILER}")
  set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_libomp_LIBRARY "omp")
  set (OpenMP_libgomp_LIBRARY "gomp")
  set (OpenMP_libiomp5_LIBRARY "iomp5")
  # and now add the OpenMP parameters to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  find_package(OpenMP REQUIRED)
  # add the flags it detects to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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

# mmmultimap (memory mapped multimap)
ExternalProject_Add(mmmultimap
  GIT_REPOSITORY "https://github.com/ekg/mmmultimap.git"
  GIT_TAG "88c734c36563048b0f3acc04dd8856f19e02b75f"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(mmmultimap SOURCE_DIR)
set(mmmultimap_INCLUDE "${SOURCE_DIR}/src")

# mmmultimap (memory mapped multimap)
ExternalProject_Add(iitii
  GIT_REPOSITORY "https://github.com/ekg/iitii.git"
  GIT_TAG "85209e07a3ee403fb6557387a7f897cd76be4406"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(iitii SOURCE_DIR)
set(iitii_INCLUDE "${SOURCE_DIR}/src")

ExternalProject_Add(mmap_allocator
  GIT_REPOSITORY "https://github.com/ekg/mmap_allocator.git"
  GIT_TAG "ed61daf094de1c2e1adbe8306287ad52da5f0264"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND "make"
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mmap_allocator SOURCE_DIR)
set(mmap_allocator_INCLUDE "${SOURCE_DIR}")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  GIT_REPOSITORY "https://github.com/SaschaWitt/ips4o.git"
  GIT_TAG "bff3ccf0bf349497f2bb10f825d160b792236367"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# BBHash perfect hasher, alltypes branch
ExternalProject_Add(bbhash
  GIT_REPOSITORY "https://github.com/vgteam/BBHash.git"
  GIT_TAG "36e4fe3eaeef762c831c49cdc01f1a3a2c7a97a4"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

ExternalProject_Add(atomicbitvector
  GIT_REPOSITORY "https://github.com/ekg/atomicbitvector.git"
  GIT_TAG "10eeaaeb3cc809c38af5255c3da42c3eeec50d1c"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicbitvector SOURCE_DIR)
set(atomicbitvector_INCLUDE "${SOURCE_DIR}/include")


set(CMAKE_BUILD_TYPE Debug)

# set up our target executable and specify its dependencies and includes
add_executable(seqwish
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  ${CMAKE_SOURCE_DIR}/src/seqindex.cpp
  ${CMAKE_SOURCE_DIR}/src/paf.cpp
  ${CMAKE_SOURCE_DIR}/src/sxs.cpp
  ${CMAKE_SOURCE_DIR}/src/cigar.cpp
  ${CMAKE_SOURCE_DIR}/src/alignments.cpp
  ${CMAKE_SOURCE_DIR}/src/pos.cpp
  ${CMAKE_SOURCE_DIR}/src/match.cpp
  ${CMAKE_SOURCE_DIR}/src/transclosure.cpp
  ${CMAKE_SOURCE_DIR}/src/links.cpp
  ${CMAKE_SOURCE_DIR}/src/compact.cpp
  ${CMAKE_SOURCE_DIR}/src/dna.cpp
  ${CMAKE_SOURCE_DIR}/src/gfa.cpp
  ${CMAKE_SOURCE_DIR}/src/vgp.cpp
  ${CMAKE_SOURCE_DIR}/src/threads.cpp
  ${CMAKE_SOURCE_DIR}/src/exists.cpp
  ${CMAKE_SOURCE_DIR}/src/mmap.cpp
  ${CMAKE_SOURCE_DIR}/src/iitii_types.cpp
  )
add_dependencies(seqwish tayweeargs)
add_dependencies(seqwish sdsl-lite)
add_dependencies(seqwish gzipreader)
add_dependencies(seqwish mmmultimap)
add_dependencies(seqwish iitii)
add_dependencies(seqwish mmap_allocator)
add_dependencies(seqwish ips4o)
add_dependencies(seqwish bbhash)
add_dependencies(seqwish atomicbitvector)
target_include_directories(seqwish PUBLIC
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${gzipreader_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${iitii_INCLUDE}"
  "${mmap_allocator_INCLUDE}"
  "${bbhash_INCLUDE}"
  "${atomicbitvector_INCLUDE}")
target_link_libraries(seqwish
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${mmap_allocator_INCLUDE}/libmmap_allocator.a"
  "-latomic"
  z)
if (BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
