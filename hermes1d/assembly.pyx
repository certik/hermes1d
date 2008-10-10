from numpy import zeros
from numpy cimport ndarray
from numpy cimport int_t as DTYPE_t


def f(ndarray[DTYPE_t] x):
    cdef ndarray y = zeros(len(x), dtype="double")
    for i in range(len(x)):
        y[i] = x[i] + 2
    return y
