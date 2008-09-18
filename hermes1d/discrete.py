class DiscreteProblem(object):

    def __init__(self):
        self.linear_form = None

    def set_num_equations(self, num):
        pass

    def set_spaces(self, space):
        self.space = space

    def set_bilinear_form(self, i, j, callback):
        self.bilinear_form = callback

    def set_linear_form(self, i, callback):
        self.linear_form = callback

    def create_matrix(self):
        from numpy import zeros
        n = self.space.ndofs()
        self.A = zeros((n, n), dtype="double")
        self.RHS = zeros((n), dtype="double")

    #@profile
    def assemble_matrix_and_rhs(self):
        count = 0
        n = len(self.space.mesh.active_elements)
        for e in self.space.mesh.iter_elements():
            count += 1
            print "done: %f%%" % (100*float(count)/n)
            shape_functions = self.space.base_functions
            n = len(shape_functions)
            for phi_i in shape_functions:
                i = phi_i.global_dof()
                for phi_j in shape_functions:
                    j = phi_j.global_dof()
                    self.A[i, j] = self.bilinear_form(phi_i, phi_j)
                if self.linear_form:
                    self.RHS[i] = self.linear_form(phi_i)
        #print self.A
        #print self.RHS
        # BC:
        penalty = 10**6
        from hermes1d import BC_DIRICHLET, BC_NEUMANN
        if self.space._bc_types[0] == BC_DIRICHLET:
            val = self.space._bc_values[0]
            n = 0
            self.A[n, n] = penalty
            self.RHS[n] = val*penalty
        if self.space._bc_types[1] == BC_DIRICHLET:
            val = self.space._bc_values[1]
            n = len(self.RHS)-1
            self.A[n, n] = penalty
            self.RHS[n] = val*penalty
        if self.space._bc_types[0] == BC_NEUMANN:
            val = self.space._bc_values[0]
            n = 0
            self.RHS[n] += val
        if self.space._bc_types[1] == BC_NEUMANN:
            val = self.space._bc_values[1]
            n = len(self.RHS)-1
            self.RHS[n] += -val

    def insert_matrix(self, mat, dof_map):
        for i in range(len(dof_map)):
            for j in range(len(dof_map)):
                self.A[dof_map[i], dof_map[j]] += mat[i, j]

    def solve_system(self, sln):
        #print self.A
        #print self.RHS
        from scipy.linalg import solve
        x = solve(self.A, self.RHS)
        #print x
        sln.set_fe_solution(self.space, x)
