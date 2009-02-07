cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(UseCython)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
set(CMAKE_CXX_FLAGS "-Wfatal-errors")

CYTHON_ADD_MODULE(assembly cassembly.cpp)
