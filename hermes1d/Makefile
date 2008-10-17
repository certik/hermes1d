#PYVER := $(shell script/config.py python_version)
#OPTFLAGS := $(shell script/config.py opt_flags)
PYVER := 2.5

CYTHON := cython --convert-range
CXX := g++
LIBS :=
CFLAGS := -Wfatal-errors -g -I/usr/include/python$(PYVER) $(OPTFLAGS) -I/usr/include/numpy -I../include
LDFLAGS := -Wfatal-errors -g

.SECONDARY:

all: assembly.so

HEADERS := cassembly.h matrix.h common.h

cassembly.o: $(HEADERS)
assembly.c: assembly.pyx assembly.pxd $(HEADERS)


%.c: %.pyx
	$(CYTHON) $<

%.o : %.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o : %.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.so: %.o cassembly.o
	$(CXX) $(LDFLAGS) -shared $+ -o $@ $(LIBS)

clean:
	rm -f *.so *.pyd *.o *.pyc assembly.c
