cmake_minimum_required(VERSION 3.10)

project(composite-measure VERSION 1.0.0)

find_package(OpenMP REQUIRED)

if(${CMAKE_CLANG})
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  set(CMAKE_CXX_COMPILER clang++)
else()
  set(CMAKE_CXX_COMPILER /usr/bin/g++)
endif()

# create library

add_subdirectory(src)

# create example executabe

add_subdirectory(example)
