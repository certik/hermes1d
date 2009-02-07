# This allows to link Cython files
# Examples:
# 1) to compile assembly.pyx to assembly.so:
#   CYTHON_ADD_MODULE(assembly)
# 2) to compile assembly.pyx and something.cpp to assembly.so:
#   CYTHON_ADD_MODULE(assembly something.cpp)

macro(CYTHON_ADD_MODULE name)
    add_custom_command(OUTPUT ${name}.cpp COMMAND cython ARGS -o ${name}.cpp ${name}.pyx DEPENDS ${name}.pyx)
    add_library(${name} SHARED ${name}.cpp ${ARGN})
    set_target_properties(${name} PROPERTIES PREFIX "")
endmacro(CYTHON_ADD_MODULE)

