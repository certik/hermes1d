class H1Space(object):

    def __init__(self, mesh, shapeset):
        self.mesh = mesh
        self.shapeset = shapeset

    def ndofs(self):
        return self.mesh.len()+1

    def get_num_dofs(self):
        return self.ndofs()

    def set_uniform_order(self, order):
        pass

    def assign_dofs(self):
        pass

    def elements(self):
        from mesh import Elem
        for e in range(self.mesh.len()):
            el = Elem(e, self.shapeset)
            el.set_dof_map((e, e+1))
            yield el
