execute_process(COMMAND python -c "import sys;print '%d.%d' % sys.version_info[:2]" OUTPUT_VARIABLE PYTHON_VERSION)
message(STATUS "XXX: ${PYTHON_VERSION}")
