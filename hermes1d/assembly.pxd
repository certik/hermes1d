cdef extern from "math.h":

    double c_sqrt "sqrt"(double x)

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void *malloc (size_t)

cdef extern from "string.h":
    void *memset(void*, int, size_t)
    void *memcpy(void*, void*, size_t)
    char *strdup(char*)

cdef extern from "arrayobject.h":

    ctypedef int intp

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags

    PyArray_EMPTY(...)
    PyArray_ZEROS(...)
    PyArray_Zeros(...)
    void* PyArray_SIZE(object)
    void* PyArray_DATA(object)

cdef extern from "Python.h":
    ctypedef void PyObject
    void Py_INCREF(PyObject *x)
    void Py_DECREF(PyObject *x)


cdef extern from "cassembly.h":

    cdef struct c_A "A":
        void set_mesh(double *mesh, int nmesh)
        void print_info()
        double *mesh
        int nmesh

    c_A *new_A "new A" ()
    void del_A "delete" (c_A *a)

cdef class A:
    cdef c_A *thisptr
