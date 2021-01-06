# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(seqwish)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Use all standard-compliant optimizations
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcx16 -g")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcx16 -g")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -mcx16 -g -fsanitize=address")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -mcx16 -g -fsanitize=address")

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
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sdsl-lite"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/args"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# gzip file reader, header only
ExternalProject_Add(gzipreader
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gzip_reader"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(gzipreader SOURCE_DIR)
set(gzipreader_INCLUDE "${SOURCE_DIR}")

# mmmulti (memory mapped multimap, multiset, and interval tree)
ExternalProject_Add(mmmulti
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmulti"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mmmulti SOURCE_DIR)
set(mmmulti_INCLUDE "${SOURCE_DIR}/src")

# iitii, like cgranges but with interpolation and memory mapping
ExternalProject_Add(iitii
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/iitii"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(iitii SOURCE_DIR)
set(iitii_INCLUDE "${SOURCE_DIR}/src")

ExternalProject_Add(mmap_allocator
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmap_allocator"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND "")
ExternalProject_Get_property(mmap_allocator INSTALL_DIR)
set(mmap_allocator_INCLUDE "${INSTALL_DIR}/include/mmap_allocator")
set(mmap_allocator_LIB "${INSTALL_DIR}/lib")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/ips4o"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# BBHash perfect hasher
ExternalProject_Add(bbhash
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/BBHash"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

# atomic bitvector class
ExternalProject_Add(atomicbitvector
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomicbitvector"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicbitvector SOURCE_DIR)
set(atomicbitvector_INCLUDE "${SOURCE_DIR}/include")

# atomic queue library
ExternalProject_Add(atomicqueue
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomic_queue"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicqueue SOURCE_DIR)
set(atomicqueue_INCLUDE "${SOURCE_DIR}")

# ska
ExternalProject_Add(ska
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/flat_hash_map"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ska SOURCE_DIR)
set(ska_INCLUDE "${SOURCE_DIR}")

# paryfor parallel_for
ExternalProject_Add(paryfor
  #GIT_REPOSITORY "https://github.com/ekg/paryfor"
  #GIT_TAG "509b28a092f732a068e2908bb9e359a8562cd32f"
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/paryfor"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(paryfor SOURCE_DIR)
set(paryfor_INCLUDE "${SOURCE_DIR}")

#add_subdirectory(deps/mmmulti/deps/mio)
ExternalProject_Add(mio
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmulti/deps/mio"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND "")
ExternalProject_Get_property(mio SOURCE_DIR)
set(mio_INCLUDE "${SOURCE_DIR}/include")

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

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
  ${CMAKE_SOURCE_DIR}/src/exists.cpp
  ${CMAKE_SOURCE_DIR}/src/time.cpp
  ${CMAKE_SOURCE_DIR}/src/mmap.cpp
  )
add_dependencies(seqwish tayweeargs)
add_dependencies(seqwish sdsl-lite)
add_dependencies(seqwish gzipreader)
add_dependencies(seqwish mmmulti)
add_dependencies(seqwish iitii)
add_dependencies(seqwish mmap_allocator)
add_dependencies(seqwish ips4o)
add_dependencies(seqwish bbhash)
add_dependencies(seqwish atomicbitvector)
add_dependencies(seqwish atomicqueue)
add_dependencies(seqwish ska)
add_dependencies(seqwish paryfor)
add_dependencies(seqwish mio)
target_include_directories(seqwish PUBLIC
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${gzipreader_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmulti_INCLUDE}"
  "${iitii_INCLUDE}"
  "${mmap_allocator_INCLUDE}"
  "${bbhash_INCLUDE}"
  "${atomicbitvector_INCLUDE}"
  "${atomicqueue_INCLUDE}"
  "${ska_INCLUDE}"
  "${paryfor_INCLUDE}"
  "${mio_INCLUDE}")
target_link_libraries(seqwish
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${mmap_allocator_LIB}/libmmap_allocator.a"
  "-latomic"
  Threads::Threads
  z)
if (BUILD_STATIC)
  #set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
endif()

install(TARGETS seqwish DESTINATION bin)
