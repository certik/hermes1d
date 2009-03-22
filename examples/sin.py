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
def F(i):
    def f1(y1, y2, t):
        return y2
    def f2(y1, y2, t):
        return -y1
    if i == 0:
        return f1
    elif i == 1:
        return f2
    raise ValueError("Wrong i (i=%d)." % (i))

# definition of the Jacobian matrix
def J(i, j):
    def f11(y1, y2, t):
        return 0
    def f12(y1, y2, t):
        return 1
    def f21(y1, y2, t):
        return -1
    def f22(y1, y2, t):
        return 0
    if i == 0 and j == 0:
        return f11
    elif i == 0 and j == 1:
        return f12
    elif i == 1 and j == 0:
        return f21
    elif i == 1 and j == 1:
        return f22
    raise ValueError("Wrong i, j (i=%d, j=%d)." % (i, j))

# assign both F and J to the discrete problem:
d.set_rhs(F, J)

# enumeration of unknowns:
d.assign_dofs()

J = d.assemble_J()
Y = zeros((J.shape[0],))
error = 1e10
i = 0
while error > 1e-1:
    F = d.assemble_F(Y)
    dY = d.solve(J, F)
    error = d.calculate_error_l2_norm(dY)
    print "it=%d, l2_norm=%e" % (i, error)
    Y += dY
    i += 1
x = Y
#print
#print J
#print F
#print x

from pylab import plot, legend, show
sln1, sln2 = d.linearize(x, 5)
x1, y1 = sln1
x2, y2 = sln2
plot(x1, y1, label="$u_1$")
plot(x2, y2, label="$u_2$")
legend()
show()
