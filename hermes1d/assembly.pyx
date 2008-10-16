from numpy import zeros, empty, float64, int32
from numpy cimport NPY_DOUBLE
from numpy cimport int_t, double_t

cimport cython

cdef inline ndarray array_i(int size, int *data):
    #cdef ndarray ary = PyArray_EMPTY(1, &size, NPY_DOUBLE, 0)
    cdef ndarray ary = empty(size, dtype=int32)
    if data != NULL: memcpy(ary.data, data, size*sizeof(int))
    return ary

cdef inline ndarray array_d(int size, double *data):
    cdef ndarray ary = empty(size, dtype=float64)
    if data != NULL: memcpy(ary.data, data, size*sizeof(double))
    return ary

cdef inline int iarray_d(ndarray a, int *size, double **data) except -1:
    if a.dtype != float64:
        raise TypeError("The array must have the dtype=float64.")
    if size!=NULL: size[0] = a.dimensions[0]
    if data!=NULL: data[0] = <double *> (a.data)

cdef extern from "cassembly.h":
    ctypedef double (*f2)(double x)
    double test2()

def test3():
    return test2()

cdef api double integ_quad_c(f2 f, double a, double b):
    from integrate import quadrature
    val, err = quadrature(func2, a, b, args=(<object>f,))
    return val

def integ_quad(f, a, b):
    from integrate import quadrature
    val, err = quadrature(func, a, b, args=(f,))
    return val

from numpy import array
def func(a, f):
    return array([f(x) for x in a])

def func2(a, f):
    cdef f2 g
    g = <f2>f
    return array([g(x) for x in a])


cdef class System:

    def __cinit__(self):
        self.thisptr = new_System()

    def __dealloc__(self):
        del_System(self.thisptr)

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

    cdef matrix2scipy(self, SparseMatrix *m):
        from scipy.sparse import coo_matrix
        Ai = array_i(m.A_len, m.Ai)
        Aj = array_i(m.A_len, m.Aj)
        Ax = array_d(m.A_len, m.Ax)
        return coo_matrix((Ax, [Ai, Aj]))

    def get_matrix_A(self):
        return self.matrix2scipy(self.thisptr.A)

    def get_matrix_B(self):
        return self.matrix2scipy(self.thisptr.B)
