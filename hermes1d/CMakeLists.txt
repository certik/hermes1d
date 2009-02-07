cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
set(CMAKE_CXX_FLAGS "-Wfatal-errors")

macro(CYTHON_ADD_MODULE name)
    add_custom_command(OUTPUT ${name}.cpp COMMAND cython ARGS -o ${name}.cpp ${name}.pyx DEPENDS ${name}.pyx)
    add_library(${name} SHARED ${name}.cpp ${ARGN})
    set_target_properties(${name} PROPERTIES PREFIX "")
endmacro(CYTHON_ADD_MODULE)

CYTHON_ADD_MODULE(assembly cassembly.cpp)
#add_custom_command(OUTPUT assembly.cpp COMMAND cython ARGS -o assembly.cpp assembly.pyx DEPENDS assembly.pyx)
#add_library(assembly SHARED cassembly.cpp assembly.cpp)
#set_target_properties(assembly PROPERTIES PREFIX "")
