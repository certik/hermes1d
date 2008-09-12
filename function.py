class Function(object):
    """
    Represents a function.

    Each function has a domain and values over this domain. However, the exact
    representation depends on the circumstances.

    Examples:

    MeshFunction() is internally represented as a linear combination of basis
    functions.

    ShapeFunction() is internally represented using the (shapeset, idx, diff)
    tuple.

    Mul() just contains the two arguments (u*v) and employs the best algorithm
    for each operation you ask it.

    The idea is that you just deal with all functions using the general
    Function interface and each function knows the best how to deal with itself
    (and others) using the most efficient algorithms.
    """

    def diff(self):
        raise NotImplementedError()

    def __call__(self, x):
        return self.f(x)

    def f(self, x):
        raise NotImplementedError()

    def domain_range(self):
        raise NotImplementedError()

    def integrate(self):
        from numpy import arange
        from scipy import trapz
        a, b = self.domain_range()
        N = 3
        x = arange(a, b, float(b-a)/N)
        y = [self.f(_x) for _x in x]
        return trapz(y, x)

    def __mul__(self, f):
        return Mul(self, f)

class Mul(Function):

    def __init__(self, a, b):
        self.args = (a, b)
        self._domain_range = a.domain_range()

    def domain_range(self):
        return self._domain_range

    def f(self, x):
        return self.args[0].f(x) * self.args[1].f(x)


class MeshFunction(Function):
    pass

class Solution(MeshFunction):

    def domain_range(self):
        return 0, 3.14159

    def f(self, x):
        idx = int(x)
        if idx < 0:
            idx = 0
        if idx > len(self.coeff) - 1:
            idx = len(self.coeff) - 1
        return self.coeff[idx]

    def set_fe_solution(self, space, x):
        self.coeff = x

class ShapeFunction(Function):

    def __init__(self, shapeset, idx, diff=0):
        self._shapeset = shapeset
        self._idx = idx
        self._diff = diff

    @property
    def idx(self):
        return self._idx

    def domain_range(self):
        return 0, 1

    def diff(self):
        return ShapeFunction(self._shapeset, self._idx, self._diff+1)

    def f(self, x):
        if self._diff == 0:
            if self._idx == 0:
                return 1-x
            if self._idx == 1:
                return x
        elif self._diff == 1:
            if self._idx == 0:
                return -1
            if self._idx == 1:
                return 1
        raise NotImplementedError()
