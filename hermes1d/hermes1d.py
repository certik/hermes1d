from numpy import empty
from quadrature import quadrature

class Node(object):
    """
    Represents a node on the mesh, given by a coordinate.
    """

    def __init__(self, x):
        self._x = x

    @property
    def x(self):
        return self._x

class Element(object):
    """
    Represents an element on the mesh, given by two nodes.
    """

    def __init__(self, x1, x2, order=1):
        self._nodes = (x1, x2)
        self._order = order
        self._dofs = [-1]*(order+1)

    @property
    def nodes(self):
        return self._nodes

    @property
    def order(self):
        return self._order

    @property
    def dofs(self):
        return self._dofs

    def assign_dofs(self, local_dofs, global_dofs):
        """
        Sets the global degrees of freedom corresponding to the local shape
        functions.

        Example 1:
        >>> e = Element(n1, n2, order=2)
        >>> e.set_dofs([0, 1, 2], [4, 331, 18])

        Example 1:
        >>> e = Element(n1, n2, order=2)
        >>> e.set_dofs([0, 1], [4, 331])
        >>> e.set_dofs([2], [18])
        """
        for l, g in zip(local_dofs, global_dofs):
            if l >= len(self._dofs):
                print l
                raise ValueError("local dof too high")
            self._dofs[l] = g

    def integrate_dphi_phi(self, i, j):
        """
        Calculates the integral of dphi*phi on the reference element
        """
        integrals_table = {
                (0, 0): -0.5,
                (1, 0): 0.5,
                (0, 1): -0.5,
                (1, 1): 0.5,
                }
        return integrals_table[(i, j)]

class Mesh(object):
    """
    Represents a finite element mesh, given by a list of nodes and then by a
    list of elements.
    """

    def __init__(self, nodes, elements):
        self._nodes = nodes
        self._elements = elements

        self._left_lift = False
        self._right_lift = False

    @property
    def nodes(self):
        return self._nodes

    @property
    def elements(self):
        return self._elements

    def set_bc(self, left=True, value=0.0):
        """
        Assign the Dirichlet bc to the mesh.

        This is needed during the assigning of dofs.

        left == True .... assign the bc to the left
        left == False ... assign to the right
        value ... the value of the function
        """
        if left:
            self._left_lift = True
            self._left_value = value
        else:
            self._right_lift = True
            self._right_value = value

    def assign_dofs(self, start_i=0):
        self._start_i = start_i
        # assign the vertex functions
        i = start_i
        if self._left_lift:
            el_list = self._elements[1:]
            self._elements[0].assign_dofs([1], [i])
        else:
            el_list = self._elements
        for e in el_list:
            e.assign_dofs([0, 1], [i, i+1])
            i += 1
        if self._right_lift:
            self._elements[-1].assign_dofs([1], [-1])
            i -= 1

        # assign bubble functions
        i += 1
        for e in self._elements:
            local_dofs = range(2, e.order+1)
            global_dofs = range(i, i+e.order-1)
            i += e.order-1
            e.assign_dofs(local_dofs, global_dofs)
        self._end_i = i
        return i

class DiscreteProblem(object):

    def __init__(self, meshes=[]):
        """
        Initializes the DiscreteProblem.

        Example:
        >>> d = DiscreteProblem(meshes=[m1, m2])
        """
        self._meshes = meshes

    def set_rhs(self, F, J):
        """
        Sets the rhs for ODE.

        Example:
        >>> e = DiscreteProblem([m1, m2])
        >>> e.set_rhs(F, J)
        """
        self._F = F
        self._J = J

    def get_mesh_number(self, global_dof_number):
        for mi, m in enumerate(self._meshes):
            if m._start_i <= global_dof_number and global_dof_number < m._end_i:
                return mi
        raise ValueError("No mesh found for global_dof=%d" % global_dof_number)

    def assign_dofs(self):
        """
        Assigns dofs for all the meshes in the problem.
        """
        i = 0
        for m in self._meshes:
            i = m.assign_dofs(start_i=i)
        self._ndofs = i
        return i

    def assemble(self):
        Y = empty((self._ndofs,))
        J = empty((self._ndofs, self._ndofs))
        for m in self._meshes:
            for e in m.elements:
                for i in range(len(e.dofs)):
                    for j in range(len(e.dofs)):
                        i_glob = e.dofs[i]
                        j_glob = e.dofs[j]
                        if i_glob == -1 or j_glob == -1:
                            continue
                        mi = self.get_mesh_number(i_glob)
                        mj = self.get_mesh_number(j_glob)
                        f = self._J(mi, mj)
                        # now f = f(y1, y2, ..., t)
                        def func(x):
                            y1 = x
                            return f(x)*phi(i, x)*phi(j, x)
                        #term_df_phi_phi = quadrature(func, -1, 1)
                        #J[i_glob, j_glob] += \
                                #    e.integrate_dphi_phi(j, i) + \
