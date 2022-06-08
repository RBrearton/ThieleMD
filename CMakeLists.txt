# minimum version for cmake
cmake_minimum_required(VERSION 3.23.1)

# project's name
project(MolecularDynamics)

# define some constant we'll use as compiler flags
SET(COMPILER_FLAGS)
add_definitions(${COMPILER_FLAGS})

# define location for produced binaries
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# NOTE: running cmake generates a makefile with explicitly named .x files, so
# every time a new .c, .cpp or .h file is added cmake should be run in order to
# generate a makefile which includes the new file.
file(GLOB SRC_FILES
  "src/*.c"
  "src/*.h"
  "src/*.cpp"
  "src/*.hpp"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
set(CMAKE_CXX_STANDARD 17)

include_directories(${X11_INCLUDE_DIR})

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)

add_executable(a.out ${SRC_FILES})

include_directories(PkgConfig::FFTW)
target_link_libraries(a.out m PkgConfig::FFTW)
