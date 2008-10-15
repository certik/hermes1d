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

    cdef struct SparseMatrix:
        int *Ai, *Aj
        double *Ax
        int A_len, A_max

    cdef struct c_System "System":
        void set_mesh(double *mesh, int nmesh)
        void print_info()
        void assemble()

        double *mesh
        int nmesh
        SparseMatrix *A, *B

    c_System *new_System "new System" ()
    void del_System "delete" (c_System *a)

cdef class System:
    cdef c_System *thisptr

    cdef matrix2numpy(self, SparseMatrix *m)
