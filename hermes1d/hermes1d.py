from math import sqrt

from quadrature import quadrature, fixed_quad

from numpy import zeros, array, arange, eye
from numpy.linalg import solve
from numpy.linalg import norm as l2_norm
from scipy.special.orthogonal import p_roots

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
    def length(self):
        return self._nodes[1].x - self._nodes[0].x

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
        def func(x):
            return self.shape_function_deriv(i, x) * self.shape_function(j, x)
        i, err = quadrature(func, -1, 1)
        return i

    def shape_function(self, idx, x):
        """
        Returns the value of the shape function "idx" at the point "x".

        "x" is in the reference domain.
        """
        if idx == 0:
            return (1-x)/2
        elif idx == 1:
            return (1+x)/2
        elif idx == 2:
            return (x**2-1)*sqrt(3./2)/2
        elif idx == 3:
            return (x**2-1)*x*sqrt(5./2)/2
        elif idx == 4:
            return (x**2-1)*(5*x**2-1)*sqrt(7./2)/8
        elif idx == 5:
            return (x**2-1)*(7*x**2-3)*sqrt(9./2)/8
        elif idx == 6:
            return (x**2-1)*(21*x**4-14*x**2+1)*sqrt(11./2)/16
        elif idx == 7:
            return (x**2-1)*(33*x**4-30*x**2+5)*sqrt(13./2)/16
        elif idx == 8:
            return (x**2-1)*(429*x**6-495*x**4+135*x**2-5)*sqrt(15./2)/128
        elif idx == 9:
            return (x**2-1)*(715*x**6-1001*x**4+385*x**2-35)*sqrt(17./2)/128
        elif idx == 10:
            return (x**2-1)*(2431*x**8-4004*x**6+2002*x**4-308*x**2+7)*sqrt(19./2)/256
        raise NotImplementedError("Such shape function is not implemented yet (i=%d)" % i)

    def shape_function_deriv(self, idx, x):
        """
        Returns the value of the derivative of the shape function "idx" at the
        point "x".

        "x" is in the reference domain.
        """
        if idx == 0:
            return -0.5
        elif idx == 1:
            return 0.5
        elif idx == 2:
            return 2*x*sqrt(3./2)/2
        elif idx == 3:
            return (3*x**2-1)*sqrt(5./2)/2
        elif idx == 4:
            return (20*x**3-12*x)*sqrt(7./2)/8
        elif idx == 5:
            return (28*x**3-20*x)*sqrt(9./2)/8
        elif idx == 6:
            return (126*x**5-140*x**3+30*x)*sqrt(11./2)/16
        elif idx == 7:
            return (198*x**5-252*x**3+70*x)*sqrt(13./2)/16
        elif idx == 8:
            return (3432*x**7-5544*x**5+2520*x**3-280*x)*sqrt(15./2)/128
        elif idx == 9:
            return (-840*x + 5544*x**3 - 10296*x**5 + 5720*x**7)*sqrt(17./2)/128
        elif idx == 10:
            return (630*x - 9240*x**3 + 36036*x**5 - 51480*x**7 + 24310*x**9)*sqrt(19./2)/256
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

    @property
    def ndofs(self):
        """
        Returns the total number of dofs.
        """
        return self._ndofs

    def assemble_J(self, Y):
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
                            nmeshes = len(self._meshes)
                            y = [0]*nmeshes
                            y1 = 0
                            y2 = 0
                            x_phys = e.ref2phys(x)
                            y.append(x_phys)
                            f_user = f(*y)
                            return f_user * \
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

    def get_sol_value(self, mesh_num, el_num, Y, x, count_lift=True):
        """
        "x" is on the *reference* element
        """
        m = self._meshes[mesh_num]
        e = m.elements[el_num]
        val = 0.
        for i, g in enumerate(e.dofs):
            if g == -1:
                if count_lift:
                    val += e.shape_function(i, x)*e.get_dirichlet_value(i)
            else:
                val += e.shape_function(i, x)*Y[g]
        #print val, e.dofs
        return val

    def assemble_F(self, Y=None):
        if Y is None:
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
                                #print "XX", e.dofs, j
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
                        y = [self.get_sol_value(_i, el_num, Y, x) for _i \
                                in range(len(self._meshes))]
                        x_phys = e.ref2phys(x)
                        y.append(x_phys)
                        return f(*y) * e.shape_function(i, x)
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
        Linearize the solution

        Y ... solution vector (all solutions)
        n ... refinement for all elements (how many refinements should be done
              to one element for the visualization purposes)

        Returns a tuple with all solutions, where each solution contains (x, y)
        points.
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

    def calculate_error_l2_norm(self, dY):
        """
        Returns the L2 norm of the vector dY.

        E.g. the square root of the sum of squares of the components of dY.
        """
        solutions = []
        norm = 0.
        for mi in range(len(self._meshes)):
            for ei in range(len(self._meshes[mi].elements)):
                e = self._meshes[mi].elements[ei]
                # change this to gauss points:
                x_vals, w = p_roots(20)
                norm_e_squared = 0.
                for i, x in enumerate(x_vals):
                    norm_e_squared += w[i] * \
                            self.get_sol_value(mi, ei, dY, x,
                                    count_lift=False)**2
                norm_e_squared *= e.jacobian
                norm += norm_e_squared
        return sqrt(norm)

    def get_initial_condition_euler(self, tol=1e-10):
        """
        Calculates the initial vector to the Newton's iteration.

        Notes:

          * This only works if the boundary conditions are all given on the
            left.  (If this is not the case, this function raises an
            exception.)
          * Nodal values are calculated using the implicit Euler method, higher
            order part is set to zero

        Todo:
          * For higher order elements one should calculate all the coefficients
            using projections
        """
        Z = zeros((len(self._meshes), len(self._meshes[0].elements)+1))
        for mi, m in enumerate(self._meshes):
            if not m._left_lift:
                raise Exception("get_initial_condition_euler() only works if all boundary conditions are given on the left.")

            Z[mi, 0] = m._left_value
        def get_F(Z, tau):
            """
            Evaluates the RHS for the vector Z and time tau.
            """
            Z0 = zeros((len(self._meshes),))
            for mi, m in enumerate(self._meshes):
                args = list(Z)+[tau]
                Z0[mi] = self._F(mi)(*args)
            return Z0
        def get_phi(Z, Zprev, tau, t):
            return Z - tau*get_F(Z, t)-tau*Zprev
        def get_J(Z, tau, t):
            mat = eye(len(self._meshes))
            for i in range(len(self._meshes)):
                for j in range(len(self._meshes)):
                    args = list(Z)+[t]
                    mat[i, j] += - tau*self._J(i, j)(*args)
            return mat

        # initial time and initial condition vector:
        tprev = self._meshes[0].elements[0].nodes[0].x
        Zprev = Z[:, 0].copy()
        Znext = Zprev[:].copy()
        for el_i in range(len(self._meshes[0].elements)):
            print "doing element:", el_i
            tau = self._meshes[0].elements[el_i].length
            print "Znext", Znext
            tnext = tprev + tau
            error = 1e10
            i = 0
            while error > tol:
                J = get_J(Zprev, tau, tprev)
                print J
                phi = get_phi(Znext, Zprev, tau, tprev)
                print phi
                stop
                dZ = solve(J, -phi)
                Znext += dZ
                error_dZ = l2_norm(dZ)
                error_phi = l2_norm(get_phi(Znext, Zprev, tau, tnext))
                print "it=%d, l2_norm_dZ=%e, l2_norm_phi=%e" %  \
                    (i, error_dZ, error_phi)
                error = max(error_dZ, error_phi)
                i += 1
            Z[:, el_i+1] = Znext[:].copy()
            Zprev = Znext[:].copy()
            tprev = tnext

        #print Z
        #stop
        #from pylab import plot, legend, show
        #plot(range(len(Z[0, :])), Z[0, :], label="$u_1$")
        #plot(range(len(Z[0, :])), Z[1, :], label="$u_2$")
        #legend()
        #show()

        #print Z
        return zeros((self._ndofs,))
