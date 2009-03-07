from numpy import zeros, array, arange
from numpy.linalg import solve
from quadrature import quadrature, fixed_quad

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
        self._lifts = [0.]*(order+1)

    @property
    def nodes(self):
        return self._nodes

    @property
    def order(self):
        return self._order

    @property
    def dofs(self):
        return self._dofs

    @property
    def jacobian(self):
        return (self._nodes[1].x - self._nodes[0].x)/2.

    def get_dirichlet_value(self, i):
        return self._lifts[i]

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

    def shape_function(self, i, x):
        """
        Returns the value of the shape function "i" at the point "x".

        "x" is in the reference domain.
        """
        if i == 0:
            return -0.5*x+0.5
        if i == 1:
            return 0.5*x+0.5
        raise NotImplementedError("Such shape function is not implemented yet (i=%d)" % i)

    def shape_function_deriv(self, i, x):
        """
        Returns the value of the shape function "i" at the point "x".

        "x" is in the reference domain.
        """
        if i == 0:
            return -0.5
        if i == 1:
            return 0.5
        raise NotImplementedError("Such shape function is not implemented yet (i=%d)" % i)

    def ref2phys(self, xi):
        """
        Converts xi to a physical domain.

        xi is in the reference domain, ref2phys returns a point in the phys.
        domain.
        """
        a = self._nodes[0].x
        b = self._nodes[1].x
        return (a+b)/2. - (a-b)/2. * xi

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
            self._elements[0]._lifts[0] = self._left_value
        else:
            el_list = self._elements
        for e in el_list:
            e.assign_dofs([0, 1], [i, i+1])
            i += 1
        if self._right_lift:
            self._elements[-1].assign_dofs([1], [-1])
            self._elements[-1]._lifts[1] = self._right_value
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

    def assemble_J(self):
        Y = zeros((self._ndofs,))
        J = zeros((self._ndofs, self._ndofs))
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
                            # x is the integration point, we need to determine
                            # the values of y1, y2, ... at this integration
                            # point.

                            # XXX: for linear problems, it doesn't matter what
                            # those values are:
                            y1 = 0
                            y2 = 0
                            x_phys = e.ref2phys(x)
                            if len(self._meshes) == 1:
                                return f(y1, x_phys) * \
                                        e.shape_function(i, x) * \
                                        e.shape_function(j, x)
                            else:
                                return f(y1, y2, x_phys) * \
                                        e.shape_function(i, x) * \
                                        e.shape_function(j, x)
                        dphi_phi = e.integrate_dphi_phi(j, i)
                        df_phi_phi, err = quadrature(func, -1, 1)
                        df_phi_phi *= e.jacobian
                        J[i_glob, j_glob] += dphi_phi - df_phi_phi
                        #print f(0, array([-1, -0.9, 0, 0.7, 0.9, 1]))
                        #stop
                        #print "X", i_glob, j_glob, i, dphi_phi, df_phi_phi
        return J

    def get_sol_value(self, mesh_num, el_num, Y, x):
        """
        "x" is on the *reference* element
        """
        m = self._meshes[mesh_num]
        e = m.elements[el_num]
        val = 0.
        for i, g in enumerate(e.dofs):
            if g == -1:
                val += e.shape_function(i, x)*e.get_dirichlet_value(i)
            else:
                val += e.shape_function(i, x)*Y[g]
        #print val, e.dofs
        return val

    def assemble_F(self):
        Y = zeros((self._ndofs,))
        F = zeros((self._ndofs,))
        for m in self._meshes:
            for el_num, e in enumerate(m.elements):
                for i in range(len(e.dofs)):
                    i_glob = e.dofs[i]
                    if i_glob == -1:
                        continue
                    mi = self.get_mesh_number(i_glob)
                    f = self._F(mi)
                    # now f = f(y1, y2, ..., t)
                    def func1(x):
                        # x is the integration point, we need to determine
                        # the values of y1, y2, ... at this integration
                        # point.

                        v = 0.
                        for j in range(len(e.dofs)):
                            g = e.dofs[j]
                            if g == -1:
                                coeff = e.get_dirichlet_value(j)
                                print "XX", e.dofs, j
                            else:
                                coeff = Y[g]
                            v += coeff*e.shape_function_deriv(j, x)
                        #print "deriv", el_num, x, v
                        v = v*e.shape_function(i, x)
                        return v
                    du_phi, err = quadrature(func1, -1, 1)
                    def func2(x):
                        # x is the integration point, we need to determine
                        # the values of y1, y2, ... at this integration
                        # point.

                        # XXX: this only works if all the meshes are the same:
                        y1 = self.get_sol_value(0, el_num, Y, x)
                        #print "y1", el_num, x, y1
                        #print x
                        x_phys = e.ref2phys(x)
                        if len(self._meshes) == 2:
                            y2 = self.get_sol_value(1, el_num, Y, x)
                        if len(self._meshes) == 1:
                            return f(y1, x_phys) * e.shape_function(i, x)
                        elif len(self._meshes) == 2:
                            return f(y1, y2, x_phys) * e.shape_function(i, x)
                    #if el_num == 0:
                    #    print "func", func2(array([-1, -0.9, -0.5, 0, 0.5, 0.9,
                    #        1]))
                    #print "f", f(0, array([-1, -0.9, -0.5, 0, 0.5, 0.9,
                    #        1]))
                    #print "func", func2(array([-1, -0.9, -0.5, 0, 0.5, 0.9,
                    #        1]))
                    #stop

                    f_phi, err = quadrature(func2, -1., 1.)
                    f_phi *= e.jacobian
                    #print "X", i_glob, el_num, i, du_phi, f_phi
                    F[i_glob] += du_phi - f_phi
        #print Y
        #print "get_sol_value"
        #print self.get_sol_value(0, 0, Y, 1)
        return F

    def solve(self, J, F):
        return solve(J, -F)

    def linearize(self, Y, n):
        """
        Y ... solution vector (all solutions)
        n ... refinement for all elements
        """
        solutions = []
        for mi in range(len(self._meshes)):
            x_list = []
            y_list = []
            for ei in range(len(self._meshes[mi].elements)):
                e = self._meshes[mi].elements[ei]
                x_vals = arange(-1, 1, 2./n)
                for x in x_vals:
                    y = self.get_sol_value(mi, ei, Y, x)
                    x_list.append(e.ref2phys(x))
                    y_list.append(y)
            solutions.append((x_list, y_list))
        return solutions
