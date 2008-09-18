from integrate import quadrature

class Mesh(object):
    """
    Mesh represents a FEM mesh with elements (Elem), but it doesn't have
    information about shapefunctions/their order (it contains helper functions
    to work with shapefunctions though).

    Use Space to represent shapefunctions.
    """

    def set_interval(self, a, b):
        self.nodes = [Node(self, a), Node(self, b)]
        self.all_elements = [Elem(self, self.nodes[0], self.nodes[1])]
        self.active_elements = self.filter_active_elements(self.all_elements)

    def iter_elements(self):
        for e in self.active_elements:
            yield e

    def filter_active_elements(self, elements):
        return [e for e in elements if e.is_active()]

    def refine_element(self, id):
        pass

    def refine_all_elements(self):
        l = []
        for e in self.active_elements:
            e1, e2, n = e.refine()
            l.append(e1)
            l.append(e2)
            self.all_elements.append(e1)
            self.all_elements.append(e2)
            self.nodes.append(n)
            e.disable()
        self.active_elements = l

    def get_num_elements(self):
        return len(self.active_elements)

    def len(self):
        return self._len

    def get_nodes_x(self):
        coords = [n.x for n in self.nodes]
        coords.sort()
        return coords

    def get_element_by_coor(self, x):
        for e in self.active_elements:
            if e.contains(x):
                return e
        # return None instead of an exception
        #raise Exception("No element contains x=%f." % x)

    def get_min_max(self):
        n1 = self.active_elements[0].nodes[0].x
        n2 = self.active_elements[-1].nodes[1].x
        return n1, n2


class Node(object):

    def __init__(self, mesh, x):
        self.mesh = mesh
        self.x = float(x)

class Elem(object):
    """
    Represents an element in a FEM mesh. Doesn't contain information about
    shapefunctions (but it has helper functions to work with them).
    """

    def __init__(self, mesh, n1, n2):
        self.mesh = mesh
        self.nodes = [n1, n2]
        self.active = True

    def set_dof_map(self, dofs):
        self.dof_map = dofs

    def is_active(self):
        return self.active

    def disable(self):
        self.active = False

    def refine(self):
        n1, n2 = self.nodes
        half = Node(self.mesh, (n1.x+n2.x)/2)
        e1 = Elem(self.mesh, self.nodes[0], half)
        e2 = Elem(self.mesh, half, self.nodes[1])
        return e1, e2, half

    def contains(self, x):
        return self.nodes[0].x <= x and x <= self.nodes[1].x

    def real2reference(self, x):
        """
        Converts "x" to local coordinates (reference element).

        -1 <= real2reference(x) <= 1
        """

        n1, n2 = self.nodes
        n1, n2 = n1.x, n2.x
        half = (n1+n2)/2
        d = n2 - n1
        return 2*(x-half)/d

    def get_jacobian(self):
        """
        Returns the Jacobian of the transformation.
        """
        a, b = self.nodes[0].x, self.nodes[1].x
        h = (b-a)/2
        J = h
        return J

    #@profile
    def integrate_function(self, f):
        """
        Integrate the function "f" on the element.
        """
        from numpy import array
        a, b = self.nodes[0].x, self.nodes[1].x
        def func(x):
            #print x
            #print array([f.f(y) for y in x])
            return array([f.f(y) for y in x])
        val, err = quadrature(func, a, b)
        #print val, a, b
        #stop
        return val
