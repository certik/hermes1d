"""
Solves the second order ODE:

y'' + y = 0
y(0)=0; y'(0)=1

So the solution is y(x) = sin(x)
"""
from hermes1d import Node, Element, Mesh, DiscreteProblem

from math import pi
from numpy import zeros

def plot_Y(Y, a, b):
    # plot the result:
    from pylab import plot, legend, show, clf, axis
    #clf()
    # linearize the solution, divide each element into 5
    sln1, sln2 = d.linearize(Y, 5)
    x1, y1 = sln1
    x2, y2 = sln2
    plot(x1, y1, label="$u_1$")
    plot(x2, y2, label="$u_2$")
    axis([a, b, -1.5, 1.5])
    #legend()
    show()

# interval end points
a = 0.
b = pi*(1./4+1)

# number of elements:
N = 20

# x values of the nodes:
x_values =[(b-a)/N * i for i in range(N+1)]

# define nodes:
nodes = [Node(x) for x in x_values]

# define elements of the 1st mesh
elements = [Element(nodes[i], nodes[i+1], order=2) for i in range(N)]
m1 = Mesh(nodes, elements)
m1.set_bc(left=True, value=0)

# define elements of the 2nd mesh
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N)]
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

# definition of the initial condition for the global Newton method:
Y = d.get_initial_condition_euler()
#plot_Y(Y, a, b)
#stop
#Y = zeros((d.ndofs,))

# Newton's iteration:
error = 1e10
i = 0
J = d.assemble_J(Y)
while error > 1e-5:
    F = d.assemble_F(Y)
    dY = d.solve(J, F)
    Y += dY
    #plot_Y(Y, a, b)
    error_dY = d.calculate_error_l2_norm(dY)
    error_F = d.calculate_error_l2_norm(d.assemble_F(Y))
    print "it=%d, l2_norm_dY=%e, l2_norm_F=%e" % (i, error_dY, error_F)
    error = max(error_dY, error_F)
    i += 1

plot_Y(Y, a, b)
