from numpy import zeros
from numpy cimport ndarray
from numpy cimport int_t, double_t

cimport cython

@cython.boundscheck(False)
def f(ndarray[double_t] x):
    cdef ndarray[double_t] y = zeros(len(x), dtype="double")
    cdef int i
    for i in range(len(x)):
        y[i] = x[i] + 2
    return y
