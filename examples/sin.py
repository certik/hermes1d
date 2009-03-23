"""
Solves the second order ODE:

y'' + y = 0
y(0)=0; y'(0)=1

So the solution is y(x) = sin(x)
"""
from hermes1d import Node, Element, Mesh, DiscreteProblem

from numpy import zeros
# interval end points
a = 0.
b = 4.

# number of elements:
N = 100

# define nodes:
nodes = [Node((b-a)/N * i) for i in range(N)]

# define elements of the 1st mesh
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N-1)]
m1 = Mesh(nodes, elements)
m1.set_bc(left=True, value=0)

# define elements of the 2nd mesh
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N-1)]
m2 = Mesh(nodes, elements)
m2.set_bc(left=True, value=1)

# definition of the ODE system:
d = DiscreteProblem(meshes=[m1, m2])

# definition of the RHS:
def F(i, Y, t):
    if i == 0:
        return Y[1]
    elif i == 1:
        return -Y[0]
    raise ValueError("Wrong i (i=%d)." % (i))

# definition of the Jacobian matrix
def DFDY(i, j, Y, t):
    if i == 0 and j == 0:
        return 0.
    elif i == 0 and j == 1:
        return 1.
    elif i == 1 and j == 0:
        return -1.
    elif i == 1 and j == 1:
        return 0.
    raise ValueError("Wrong i, j (i=%d, j=%d)." % (i, j))

# assign both F and J to the discrete problem:
d.define_ode(F, DFDY)

# enumeration of unknowns:
d.assign_dofs()

# definition of the initial condition for the Newton method:
Y = d.get_initial_condition_euler()
#Y = zeros((d.ndofs,))

# Newton's iteration:
error = 1e10
i = 0
while error > 1e-5:
    F = d.assemble_F(Y)
    J = d.assemble_J(Y)
    dY = d.solve(J, F)
    error_dY = d.calculate_error_l2_norm(dY)
    Y += dY
    error_F = d.calculate_error_l2_norm(d.assemble_F(Y))
    print "it=%d, l2_norm_dY=%e, l2_norm_F=%e" % (i, error_dY, error_F)
    error = max(error_dY, error_F)
    i += 1

# plot the result:
from pylab import plot, legend, show
# linearize the solution, divide each element into 5
sln1, sln2 = d.linearize(Y, 5)
x1, y1 = sln1
x2, y2 = sln2
plot(x1, y1, label="$u_1$")
plot(x2, y2, label="$u_2$")
legend()
show()
