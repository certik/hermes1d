class Mesh(object):

    def set_interval(self, a, b):
        self.domain_range = (a, b)
        self._len = 2

    def refine_element(self, id):
        pass

    def refine_all_elements(self):
        self._len *= 2

    def len(self):
        return self._len

    def get_nodes_x(self):
        from numpy import arange
        a, b = self.domain_range
        return arange(a, b, float(b-a)/(self._len+1))

class Elem(object):

    def __init__(self, id, shapeset):
        self.shapeset = shapeset

    def set_dof_map(self, dofs):
        self.dof_map = dofs

    def shape_functions(self):
        from function import ShapeFunction
        for idx in range(len(self.dof_map)):
            yield ShapeFunction(self.shapeset, idx)
