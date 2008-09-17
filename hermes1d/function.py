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

    def is_zero(self, tol=1e-6):
        """
        Compares "self" with a zero function.

        Currently this is only done by comparing node values.
        """
        for n in self.mesh.nodes:
            if self.f(n.x) > tol:
                return False
        return True

    def __mul__(self, f):
        f = convert(f, self.mesh)
        return Mul(self, f)

    def __rmul__(self, f):
        f = convert(f, self.mesh)
        return Mul(f, self)

    def __pow__(self, a):
        a = convert(a, self.mesh)
        return Pow(self, a)

    def __neg__(self):
        return (-1) * self

    def __add__(self, a):
        a = convert(a, self.mesh)
        return Add(self, a)

    def __radd__(self, a):
        a = convert(a, self.mesh)
        return Add(a, self)

    def __sub__(self, a):
        a = convert(a, self.mesh)
        return Add(self, -a)

    def __rsub__(self, a):
        a = convert(a, self.mesh)
        return Add(a, -self)

    def __div__(self, a):
        a = convert(a, self.mesh)
        return self * (a ** (-1))

def convert(a, mesh=None):
    if isinstance(a, Function):
        return a
    elif isinstance(a, (int, long, float)):
        return ConstantFunction(mesh, a)
    raise NotImplementedError("Don't know how to convert to Function.")

class Pow(Function):

    def __init__(self, a, b):
        self.args = (a, b)
        self.mesh = a.mesh

    def domain_elements(self):
        return self.args[0].domain_elements()

    def f(self, x):
        a, b = self.args
        return a.f(x) ** b.f(x)

class Add(Function):

    def __init__(self, a, b):
        self.args = (a, b)
        self.mesh = a.mesh

    def f(self, x):
        a, b = self.args
        return a.f(x) + b.f(x)

class Mul(Function):

    def __init__(self, a, b):
        self.args = (a, b)
        self.mesh = a.mesh

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

class ConstantFunction(Function):
    """
    Represents a constant function.
    """

    def __init__(self, mesh, c):
        self._c = c
        self.mesh = mesh

    def f(self, x):
        return self._c

class CustomFunction(Function):
    """
    Represents a constant function.
    """

    def __init__(self, F, space):
        self._F = F
        self.space = space
        self.mesh = space.mesh

    def domain_elements(self):
        return self.mesh.active_elements

    def f(self, x):
        return self._F(x)

class LinearFunction(Function):
    """
    Represents an "x" function, i.e. y(x) = x.
    """

    def __init__(self, mesh):
        self.mesh = mesh

    def f(self, x):
        return x

    def domain_elements(self):
        return self.mesh.active_elements


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
        from numpy import arange
        x0 = []
        y0 = []
        for e in self.els:
            idx = self.els[e]
            if idx in [0, 1]:
                steps = 3
            else:
                steps = idx*10
            a = e.nodes[0].x
            b = e.nodes[1].x
            x = arange(a, b, (b-a)/steps)
            #print x[1]
            #print self.mesh.get_element_by_coor(x[1]), e
            y = [self.values(e.real2reference(xx), idx, 1) for xx in x]
            x0.extend(x)
            y0.extend(y)
        return x0, y0

    def eval_deriv(self, x, order=1):
        e = self.mesh.get_element_by_coor(x)
        if e is not None and e in self.els:
            if order == 0:
                J = 1.
            elif order == 1:
                J = 1./e.get_jacobian()
            elif order == 2:
                # XXX: verify:
                J = 1./e.get_jacobian()**2
            else:
                raise NotImplementedError()
            #print J*self.values(e.real2reference(x), self.els[e], order)
            return J*self.values(e.real2reference(x), self.els[e], order)
        else:
            return 0.
