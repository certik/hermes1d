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

    def get_xy(self, steps=5):
        def interp(a, b, steps):
            x = [a]
            dx = float(b-a)/steps
            while x[-1] + dx < b:
                x.append(x[-1]+dx)
            return x

        x = []
        for e in self.mesh.iter_elements():
            x.extend(interp(e.nodes[0].x, e.nodes[1].x, steps))
        x.sort()
        y = [self.f(c) for c in x]
        return x, y


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
        self.mesh = f.mesh

    def f(self, x):
        return self._f.eval_deriv(x)

    def domain_elements(self):
        return self._f.domain_elements()

class MeshFunction(Function):
    """
    Represent any function that is defined as a linear combination of
    BaseFunctions.
    """

    #def f(self, x):
    #    e = self.mesh.get_element_x(x)
    #    return self.coeff[e.id]

class Solution(MeshFunction):

    def domain_range(self):
        return self.mesh.domain_range

    def set_fe_solution(self, space, x):
        self.space = space
        self.mesh = space.mesh
        self.coeff = x

    def f(self, x):
        val = 0.
        for c, b in zip(self.coeff, self.space.base_functions):
            val += c*b.f(x)
        return val

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
        self.els = {}

    def global_dof(self):
        return self.dof

    def contains(self, el):
        return el in self.els

    def domain_elements(self):
        return [e for e, idx in self.els.iteritems()]

    def element_idx(self, el):
        return self.els[el]

    def add_element(self, e, idx):
        """
        Defines the BaseFunction on the element "e".

        e ... Elem instance
        idx ... integer, idx of the shape function from the shapeset

        """
        self.els[e] = idx

    def values(self, x, idx, diff=0):
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

    def f(self, x):
        e = self.mesh.get_element_by_coor(x)
        if e in self.els:
            return self.values(e.real2reference(x), self.els[e], 0)
        else:
            return 0.

    def eval_deriv(self, x):
        e = self.mesh.get_element_by_coor(x)
        if e in self.els:
            return self.values(e.real2reference(x), self.els[e], 1)
        else:
            return 0.
