class Function(object):
    """
    Represents a function.

    Each function has a domain and values over this domain. However, the exact
    representation depends on the circumstances.

    Examples:

    MeshFunction() is internally represented as a linear combination of basis
    functions.

    BaseFunction() is internally represented using the (shapeset, idx, diff)
    tuple for all elements over which it spans and by some particular mesh.

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
    """
    Represent any function that is defined as a linear combination of
    BaseFunctions.
    """

    def get_xy(self, n=3):
        """
        Returns linearized xy values.

        n ... the number of points between mesh elements.
        """
        return self.mesh.get_nodes_x(), self.coeff

    def f(self, x):
        e = self.mesh.get_element_x(x)
        return self.coeff[e.id]

class Solution(MeshFunction):

    def domain_range(self):
        return self.mesh.domain_range

    def set_fe_solution(self, space, x):
        self.space = space
        self.mesh = space.mesh
        self.coeff = x

class xShapeFunction(Function):
    """
    Deprecated.
    """

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

class BaseFunction(Function):
    """
    Represents one base function.

    It is internally represented using the (shapeset, idx, diff) tuple for all
    elements over which it spans and by some particular mesh.
    """

    def __init__(self, mesh, shapeset):
        self.mesh = mesh
        self.shapeset = shapeset
        self.els = []

    def get_xy(self, steps=2):
        d = {}
        for e, idx in self.els:
            d[e] = idx
        x = []
        y = []
        for e in self.mesh.iter_elements():
            x.extend([e.nodes[0].x, e.nodes[1].x])
            if e in d:
                idx = d[e]
                if idx == 0:
                    y.extend([1, 0])
                elif idx == 1:
                    y.extend([0, 1])
                else:
                    raise NotImplementedError()
            else:
                y.extend([0, 0])

        return x, y

    def add_element(self, e, idx):
        """
        Defines the BaseFunction on the element "e".

        e ... Elem instance
        idx ... integer, idx of the shape function from the shapeset

        """
        self.els.append((e, idx))
