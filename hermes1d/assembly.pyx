from numpy import zeros
from numpy cimport ndarray
from numpy cimport int_t as DTYPE_t

cimport cython

@cython.boundscheck(False)
def f(ndarray[DTYPE_t] x):
    cdef ndarray[DTYPE_t] y = zeros(len(x), dtype="double")
    cdef int i
    for i in range(len(x)):
        y[i] = x[i] + 2
    return y
