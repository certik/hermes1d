class DiscreteProblem:

    def set_num_equations(self, num):
        pass

    def set_spaces(self, space):
        pass

    def set_bilinear_form(self, i, j, callback):
        pass

    def set_linear_form(self, i, callback):
        pass

    def create_matrix(self):
        pass

    def assemble_matrix_and_rhs(self):
        pass

    def solve_system(self, sln):
        pass
