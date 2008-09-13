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
        return Derivative(self)

    def __call__(self, x):
        return self.f(x)

    def f(self, x):
        raise NotImplementedError()

    def domain_range(self):
        raise NotImplementedError()

    def integrate(self):
        from numpy import arange
        from scipy import trapz
        domain = self.domain_elements()
        integral = 0
        for e in domain:
            N = 3
            a, b = e.nodes[0].x, e.nodes[1].x
            x = arange(a, b, float(b-a)/N)
            y = [self.f(_x) for _x in x]
            integral += trapz(y, x)
        return integral

    def __mul__(self, f):
        return Mul(self, f)

class Mul(Function):

    def __init__(self, a, b):
        self.args = (a, b)

    def domain_elements(self):
        return self.args[0].domain_elements()

    def f(self, x):
        return self.args[0].f(x) * self.args[1].f(x)

class Derivative(Function):
    """
    Represents a derivative.
    """

    def __init__(self, f):
        self._f = f

    def f(self, x):
        return self._f.eval_deriv(x)

    def domain_elements(self):
        return self._f.domain_elements()

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

class BaseFunction(Function):
    """
    Represents one base function.

    It is internally represented using the (shapeset, idx, diff) tuple for all
    elements over which it spans and by some particular mesh.
    """

    def __init__(self, mesh, shapeset, dof):
        self.dof = dof
        self.mesh = mesh
        self.shapeset = shapeset
        self.els = []

    def global_dof(self):
        return self.dof

    def contains(self, el):
        d = {}
        for e, idx in self.els:
            d[e] = idx
        return el in d

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

    def domain_elements(self):
        return [e for e, idx in self.els]

    def element_idx(self, el):
        d = {}
        for e, idx in self.els:
            d[e] = idx
        assert el in d
        return d[el]

    def add_element(self, e, idx):
        """
        Defines the BaseFunction on the element "e".

        e ... Elem instance
        idx ... integer, idx of the shape function from the shapeset

        """
        self.els.append((e, idx))

    def fx(self, x, idx, diff=0):
        if diff == 0:
            if idx == 0:
                return 1-x
            if idx == 1:
                return x
        elif diff == 1:
            if idx == 0:
                return -1
            if idx == 1:
                return 1
        raise NotImplementedError()

    def fd(self, x, diff=0):
        d = {}
        for e, idx in self.els:
            d[e] = idx
        e = self.mesh.get_element_by_coor(x)
        if e in d:
            idx = d[e]
            return self.fx(x, idx, diff)
        else:
            return 0.

    def f(self, x):
        return self.fd(x, 0)

    def eval_deriv(self, x):
        return self.fd(x, 1)
