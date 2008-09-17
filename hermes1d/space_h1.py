BC_NEUMANN = 0
BC_DIRICHLET = 1

class H1Space(object):
    """
    Represents all the basis (shape) functions + their order.
    It is just a set of BaseFunctions.

    It is constructed from mesh+shapeset. This class contains all information
    needed to evaluate basis functions at any point.

    Each base (shape) function is represented using BaseFunction, which is a
    subclass of Function, so you can manipulate with it as with any other
    Function.
    """

    def __init__(self, mesh, shapeset):
        self.mesh = mesh
        self.shapeset = shapeset
        self.base_functions = None
        self.orders = None
        self._bc_types = [None, None]
        self._bc_values = [None, None]

    def ndofs(self):
        if self.base_functions is None:
            raise Exception("Call build() first.")
        return len(self.base_functions)

    def get_num_dofs(self):
        return self.ndofs()

    def set_uniform_order(self, order):
        self.orders = [order]*self.mesh.get_num_elements()

    def set_bc_types(self, a, b):
        self._bc_types = [a, b]

    def set_bc_values(self, a, b):
        self._bc_values = [a, b]

    def build(self):
        """
        Constructs all the BaseFunctions using self.mesh, self.shapeset and
        given orders on each element.
        """
        from hermes1d import BaseFunction
        b = []
        dof = 0
        f = BaseFunction(self.mesh, self.shapeset, dof)
        for i, e in enumerate(self.mesh.iter_elements()):
            order = self.orders[i]
            assert order >= 1
            f.add_element(e, 0)
            b.append(f)
            dof += 1
            f = BaseFunction(self.mesh, self.shapeset, dof)
            f.add_element(e, 1)
            for ii in range(2, order+1):
                dof += 1
                bubble = BaseFunction(self.mesh, self.shapeset, dof)
                bubble.add_element(e, ii)
                b.append(bubble)
        #del b[0]
        b.append(f)
        self.base_functions = b

    def get_base_function(self, idx):
        return self.base_functions[idx]

    def assign_dofs(self):
        self.build()

    def xelements(self):
        from mesh import Elem
        for e in range(self.mesh.len()):
            el = Elem(e, self.shapeset)
            el.set_dof_map((e, e+1))
            yield el

    def shape_functions(self, e):
        l = []
        for f in self.base_functions:
            if f.contains(e):
                l.append(f)
        return l

    def dof_map(self, e):
        l = self.shape_functions(e)
        dofmap = [b.global_dof() for b in l]
        return dofmap

