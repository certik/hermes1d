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

    def ndofs(self):
        if self.base_functions is None:
            raise Exception("Call build() first.")
        return len(self.base_functions)

    def get_num_dofs(self):
        return self.ndofs()

    def set_uniform_order(self, order):
        self.orders = [order]*self.mesh.get_num_elements()

    def build(self):
        """
        Constructs all the BaseFunctions using self.mesh, self.shapeset and
        given orders on each element.
        """
        from hermes1d import BaseFunction
        b = []
        f = BaseFunction(self.mesh, self.shapeset)
        for e in self.mesh.iter_elements():
            f.add_element(e, 0)
            b.append(f)
            f = BaseFunction(self.mesh, self.shapeset)
            f.add_element(e, 1)
        b.append(f)
        self.base_functions = b

    def get_base_function(self, idx):
        return self.base_functions[idx]

    def assign_dofs(self):
        self.build()

    def elements(self):
        from mesh import Elem
        for e in range(self.mesh.len()):
            el = Elem(e, self.shapeset)
            el.set_dof_map((e, e+1))
            yield el
