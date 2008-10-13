#PYVER := $(shell script/config.py python_version)
#OPTFLAGS := $(shell script/config.py opt_flags)
PYVER := 2.5

CYTHON := cython --convert-range
CXX := g++
LIBS :=
CFLAGS := -Wfatal-errors -I/usr/include/python$(PYVER) $(OPTFLAGS) -I/usr/include/numpy -I../include
LDFLAGS := -Wfatal-errors

.SECONDARY:

all: assembly.so

%.c: %.pyx
	$(CYTHON) $<

%.o : %.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o : %.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

cassembly.o: cassembly.h
assembly.c: assembly.pyx assembly.pxd

%.so: %.o cassembly.o
	$(CXX) $(LDFLAGS) -shared $+ -o $@ $(LIBS)

clean:
	rm -f *.so *.pyd *.o *.pyc assembly.c
