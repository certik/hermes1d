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
        domain = self.domain_elements()
        integral = 0
        for e in domain:
            integral += e.integrate_function(self)
        return integral

    def __mul__(self, f):
        return Mul(self, f)

class Mul(Function):

    def __init__(self, a, b):
        self.args = (a, b)

    def domain_elements(self):
        s = None
        for x in self.args:
            if s is None:
                s = set(x.domain_elements())
            else:
                s = s.intersection(set(x.domain_elements()))
        s = list(s)
        return s

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

    def eval_deriv(self, x, order=1):
        return self._f.eval_deriv(x, 2)

    def domain_elements(self):
        return self._f.domain_elements()

class MeshFunction(Function):
    """
    Represent any function that is defined as a linear combination of
    BaseFunctions.
    """
    pass

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

    def eval_deriv(self, x, order=1):
        val = 0.
        for c, b in zip(self.coeff, self.space.base_functions):
            val += c*b.eval_deriv(x, order)
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
        return self.shapeset.get_value_reference(x, idx, diff)

    def f(self, x):
        # order=0 means the function value:
        return self.eval_deriv(x, order=0)

    def get_xy(self, steps=5):
        def interp(a, b, steps):
            x = [a]
            dx = float(b-a)/steps
            for i in range(steps):
                x.append(x[-1]+dx)
            return x

        x = []
        for e in self.mesh.iter_elements():
            if e in self.els:
                idx = self.els[e]
                if idx in [0, 1]:
                    steps = 1
                else:
                    steps = idx*10
                x.extend(interp(e.nodes[0].x, e.nodes[1].x, steps))
            else:
                x.extend((e.nodes[0].x, e.nodes[1].x))
        x.sort()
        y = [self.f(c) for c in x]
        return x, y

    def eval_deriv(self, x, order=1):
        e = self.mesh.get_element_by_coor(x)
        if e is not None and e in self.els:
            J = 1.
            if order == 1:
                # XXX: make this more general:
                # this is needed because derivatives need to be transformed.
                a, b = e.nodes[0].x, e.nodes[1].x
                h = b-a
                J = 1/h
                stop
            return J*self.values(e.real2reference(x), self.els[e], order)
        else:
            return 0.
