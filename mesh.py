class Mesh(object):

    def set_interval(self, a, b):
        pass

    def refine_element(self, id):
        pass

    def len(self):
        return 5

class Elem(object):

    def __init__(self, id, shapeset):
        self.shapeset = shapeset

    def set_dof_map(self, dofs):
        self.dof_map = dofs

    def shape_functions(self):
        from function import ShapeFunction
        for idx in range(len(self.dof_map)):
            yield ShapeFunction(self.shapeset, idx)
