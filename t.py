#! /usr/bin/env python

from hermes1d import finalize, Mesh, H1Shapeset, H1Space, DiscreteProblem, \
        Solution, ScalarView, int_grad_u_grad_v, int_u_v, BaseView, MeshView, \
        MatrixView, BC_DIRICHLET, BC_NEUMANN, LinearFunction, BaseFunction

#def bilinear_form(u, v):
#    return -int_grad_u_grad_v(u, v)+int_u_v(u, v)

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
mesh.refine_all_elements()
mesh.refine_all_elements()
shapeset = H1Shapeset()

print "Creating base functions"
space = H1Space(mesh, shapeset)
space.set_uniform_order(1)
space.set_bc_types(BC_DIRICHLET, BC_DIRICHLET)
space.set_bc_values(0, 0)
space.build()
#v = ScalarView()
#b = BaseFunction(mesh, shapeset, 1)
#e = mesh.all_elements[2]
#b.add_element(e, 0)
#e = mesh.all_elements[3]
#b.add_element(e, 1)
#v.show(LinearFunction(mesh))
#v.show(b)
#bview = BaseView("Base", 0, 0, 400, 400)
#bview.show(space)
#finalize()

print "Setting up forms"
dp = DiscreteProblem()
dp.set_num_equations(1)
dp.set_spaces(space)
dp.set_bilinear_form(0, 0, bilinear_form)
dp.set_linear_form(0, linear_form)

from hermes1d.assembly import f, A
from numpy import array
s = array([1, 4.3])
#print f(s)

a = A()
a.set_mesh(s)
a.print_info()
print a.get_mesh()


#sln = Solution()
#dp.create_matrix()
#print "Assembling"
#dp.assemble_matrix_and_rhs()
#print "Solving"
#dp.solve_system(sln)

#print "Plotting"
# visualize the solution
#w = 400
#h = 300
#mview = MeshView("Mesh", 0, 0, w, h)
#mview.show(mesh)
#bview = BaseView("Base", w, 0, w, h)
#bview.show(space)
#view = ScalarView("Solution", 2*w, 0, w, h)
#view.show(sln)
##view = ScalarView("Solution", 2*w, h, w, h)
#x = LinearFunction(mesh)
#sln_exact = (1-x**2)
#error = (sln-sln_exact)**2
#view.show(error)
#print error.is_zero(1e-9)
#mat1 = MatrixView("Matrix A (%dx%d)" % dp.A.shape, 0, h, w, h)
#mat1.show(dp)

#view = ScalarView("sol'", w, h, w, h)
#view.show(sln.diff())
#view = ScalarView("sol''", 2*w, h, w, h)
#view.show(sln.diff().diff())

#finalize()
