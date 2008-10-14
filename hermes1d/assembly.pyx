from numpy import zeros, float64
from numpy cimport NPY_DOUBLE
from numpy cimport int_t, double_t

cimport cython

cdef inline ndarray array_d(int size, double *data):
    #cdef ndarray ary = PyArray_EMPTY(1, &size, NPY_DOUBLE, 0)
    cdef ndarray ary = zeros(size, dtype=float64)
    if data != NULL: memcpy(ary.data, data, size*sizeof(double))
    return ary

cdef inline int iarray_d(ndarray a, int *size, double **data) except -1:
    if a.dtype != float64:
        raise TypeError("The array must have the dtype=float64.")
    if size!=NULL: size[0] = a.dimensions[0]
    if data!=NULL: data[0] = <double *> (a.data)


@cython.boundscheck(False)
def f(ndarray[double_t] x):
    cdef ndarray[double_t] y = zeros(len(x), dtype="double")
    cdef int i
    for i in range(len(x)):
        y[i] = x[i] + 2
    return y


cdef class A:

    def __cinit__(self):
        self.thisptr = new_A()

    def __dealloc__(self):
        del_A(self.thisptr)

    def set_mesh(self, ndarray a):
        cdef double *b
        cdef int size
        iarray_d(a, &size, &b)
        self.thisptr.set_mesh(b, size)

    def get_mesh(self):
        return array_d(self.thisptr.nmesh, self.thisptr.mesh)

    def print_info(self):
        self.thisptr.print_info()

    def assemble(self):
        self.thisptr.assemble()
