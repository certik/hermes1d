cmake_minimum_required(VERSION 2.6)
set(_CURRENT_VERSION 2.5)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PythonInterp)
find_package(PythonLibs)

include_directories(${PYTHON_INCLUDE_PATH})

add_library(assembly cassembly.cpp)
