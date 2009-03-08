project(hermes1d)

cmake_minimum_required(VERSION 2.6)

if(EXISTS CMake.vars)
    include(CMake.vars)
endif(EXISTS CMake.vars)

add_subdirectory(hermes1d)
