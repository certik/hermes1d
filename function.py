class Function(object):
    """
    Represents a function.

    Each function has a domain and values over this domain. However, the exact
    representation depends on circumstances.
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

    def __mul__(self, x):
        return self

class Solution(Function):

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
