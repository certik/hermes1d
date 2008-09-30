#! /usr/bin/env python

from hermes1d import finalize, Mesh, H1Shapeset, H1Space, DiscreteProblem, \
        Solution, ScalarView, int_grad_u_grad_v, int_u_v, BaseView, MeshView, \
        MatrixView, BC_DIRICHLET, BC_NEUMANN, LinearFunction

def test_ifem1():
    def bilinear_form(u, v):
        return (u.diff()*v.diff()).integrate()

    def linear_form(v):
        return v.integrate()

    print "Creating mesh"
    mesh = Mesh()
    mesh.set_interval(-1, 1)
    #mesh.refine_element(0)
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    shapeset = H1Shapeset()

    print "Creating base functions"
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(1)
    space.set_bc_types(BC_DIRICHLET, BC_DIRICHLET)
    space.set_bc_values(0, 0)
    space.build()
    #bview = BaseView("Base", 0, 0, 400, 400)
    #bview.show(space)
    #finalize()

    print "Setting up forms"
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_bilinear_form(0, 0, bilinear_form)
    dp.set_linear_form(0, linear_form)

    sln = Solution()
    dp.create_matrix()
    print "Assembling"
    dp.assemble_matrix_and_rhs()
    print "Solving"
    dp.solve_system(sln)
    print "comparing"

    x = LinearFunction(mesh)
    sln_exact = (1-x**2)/2
    error = (sln-sln_exact)**2
    assert error.is_zero(1e-9)

def test_ifem2():
    def bilinear_form(u, v):
        return (u.diff()*v.diff()).integrate()

    def linear_form(v):
        return 2*v.integrate()

    print "Creating mesh"
    mesh = Mesh()
    mesh.set_interval(-1, 1)
    #mesh.refine_element(0)
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    shapeset = H1Shapeset()

    print "Creating base functions"
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(1)
    space.set_bc_types(BC_DIRICHLET, BC_DIRICHLET)
    space.set_bc_values(0, 0)
    space.build()
    #bview = BaseView("Base", 0, 0, 400, 400)
    #bview.show(space)
    #finalize()

    print "Setting up forms"
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_bilinear_form(0, 0, bilinear_form)
    dp.set_linear_form(0, linear_form)

    sln = Solution()
    dp.create_matrix()
    print "Assembling"
    dp.assemble_matrix_and_rhs()
    print "Solving"
    dp.solve_system(sln)
    print "comparing"

    x = LinearFunction(mesh)
    sln_exact = (1-x**2)
    error = (sln-sln_exact)**2
    assert error.is_zero(1e-9)

def test_ifem3():
    def bilinear_form(u, v):
        return -(u.diff()*v.diff()).integrate()

    print "Creating mesh"
    mesh = Mesh()
    mesh.set_interval(-1, 1)
    #mesh.refine_element(0)
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    shapeset = H1Shapeset()

    print "Creating base functions"
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(1)
    space.set_bc_types(BC_DIRICHLET, BC_NEUMANN)
    space.set_bc_values(0, 1)
    space.build()

    print "Setting up forms"
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_bilinear_form(0, 0, bilinear_form)

    sln = Solution()
    dp.create_matrix()
    print "Assembling"
    dp.assemble_matrix_and_rhs()
    print "Solving"
    dp.solve_system(sln)
    print "comparing"

    x = LinearFunction(mesh)
    sln_exact = 1+x
    error = (sln-sln_exact)**2
    assert error.is_zero(1e-9)

def test_ifem4():
    def bilinear_form(u, v):
        return -(u.diff()*v.diff()).integrate()

    print "Creating mesh"
    mesh = Mesh()
    mesh.set_interval(-1, 1)
    #mesh.refine_element(0)
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    mesh.refine_all_elements()
    shapeset = H1Shapeset()

    print "Creating base functions"
    space = H1Space(mesh, shapeset)
    space.set_uniform_order(1)
    space.set_bc_types(BC_NEUMANN, BC_DIRICHLET)
    space.set_bc_values(1, 0)
    space.build()

    print "Setting up forms"
    dp = DiscreteProblem()
    dp.set_num_equations(1)
    dp.set_spaces(space)
    dp.set_bilinear_form(0, 0, bilinear_form)

    sln = Solution()
    dp.create_matrix()
    print "Assembling"
    dp.assemble_matrix_and_rhs()
    print "Solving"
    dp.solve_system(sln)
    print "comparing"

    x = LinearFunction(mesh)
    sln_exact = x-1
    error = (sln-sln_exact)**2
    assert error.is_zero(1e-9)
    sln_wrong = x-1.1
    error = (sln-sln_wrong)**2
    assert not error.is_zero(1e-9)
