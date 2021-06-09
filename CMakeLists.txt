cmake_minimum_required(VERSION 3.10)

project(composite-measure VERSION 1.0.0)

find_package(OpenMP REQUIRED)

find_package(PythonInterp REQUIRED)

find_package(PythonLibs REQUIRED)

include(GNUInstallDirs)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(SKBUILD)

find_package(PythonExtensions REQUIRED)

find_package(Cython REQUIRED)

find_package(NumPy REQUIRED)

endif()
message(STATUS ${PythonLibs_LIBRARIES})

if(${CMAKE_CLANG})
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  set(CMAKE_CXX_COMPILER clang++)
else()
  set(CMAKE_CXX_COMPILER /usr/bin/g++)
endif()

# create library

add_subdirectory(src)

# create example executable

add_subdirectory(example)

# create cython wraper

if(SKBUILD)
add_subdirectory(cython-wraper)
endif()
