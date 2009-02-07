cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PythonLibs)

include_directories(${PYTHON_INCLUDE_PATH})
set(CMAKE_CXX_FLAGS "-Wfatal-errors")

add_custom_command(OUTPUT assembly.cpp COMMAND cython ARGS -o assembly.cpp assembly.pyx DEPENDS assembly.pyx)
add_library(assembly cassembly.cpp assembly.cpp)
set_target_properties(assembly PROPERTIES PREFIX "")
