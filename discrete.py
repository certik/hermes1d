class DiscreteProblem:

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

    def assemble_matrix_and_rhs(self):
        from numpy import zeros
        for e in self.space.elements():
            mat = zeros((2, 2), dtype="double")
            for phi_i in e.shape_functions():
                for phi_j in e.shape_functions():
                    mat[phi_i.idx, phi_j.idx] = self.bilinear_form(phi_i,
                            phi_j)
            self.insert_matrix(mat, e.dof_map)
            #self.insert_vec(mat, e.dof_map)
        # BC:
        l = 0
        r = len(self.RHS)-1
        penalty = 10**6
        val1 = 0.0001
        val2 = 1.0
        self.A[l, l] = penalty
        self.RHS[l] = val1*penalty
        self.A[r, r] = penalty
        self.RHS[r] = val2*penalty

    def insert_matrix(self, mat, dof_map):
        for i in range(len(dof_map)):
            for j in range(len(dof_map)):
                self.A[dof_map[i], dof_map[j]] = mat[i, j]

    def solve_system(self, sln):
        print self.A
        print self.RHS
        from scipy.linalg import solve
        x = solve(self.A, self.RHS)
        print x
        sln.set_fe_solution(self.space, x)
