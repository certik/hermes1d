class H1Space(object):

    def __init__(self, mesh, shapeset):
        self.mesh = mesh
        self.shapeset = shapeset

    def set_uniform_order(self, order):
        pass

    def assign_dofs(self):
        pass

    def elements(self):
        from mesh import Elem
        for e in range(self.mesh.len()):
            el = Elem(e, self.shapeset)
            if e == 0:
                el.set_dofs((0, 1))
            elif e == 1:
                el.set_dofs((1, 2))
            else:
                raise NotImplementedError()
            yield el
