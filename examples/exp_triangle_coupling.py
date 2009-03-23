"""
Solves the system of ODE, that are coupled by a triangular matrix.
"""
from hermes1d import Node, Element, Mesh, DiscreteProblem

from numpy import zeros
a = 0.
b = 1.
N = 10
nodes = [Node((b-a)/N * i) for i in range(N)]
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N-1)]
m1 = Mesh(nodes, elements)
m1.set_bc(left=True, value=1)
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N-1)]
m2 = Mesh(nodes, elements)
m2.set_bc(left=True, value=2)

d = DiscreteProblem(meshes=[m1, m2])
k = 1.0
def F(i, Y, t):
    if i == 0:
        return -Y[0]+Y[1]
    elif i == 1:
        return Y[1]
    raise ValueError("Wrong i (i=%d)." % (i))
def DFDY(i, j, Y, t):
    if i == 0 and j == 0:
        return -1
    elif i == 0 and j == 1:
        return 1
    elif i == 1 and j == 0:
        return 0
    elif i == 1 and j == 1:
        return 1
    raise ValueError("Wrong i, j (i=%d, j=%d)." % (i, j))
d.define_ode(F, DFDY)
d.assign_dofs()
Y = zeros(d.ndofs)
J = d.assemble_J(Y)
error = 1e10
i = 0
while error > 1e-10:
    F = d.assemble_F(Y)
    dY = d.solve(J, F)
    error = d.calculate_error_l2_norm(dY)
    print "it=%d, l2_norm=%e" % (i, error)
    Y += dY
    i += 1
x = Y

from pylab import plot, legend, show
sln1, sln2 = d.linearize(x, 5)
x1, y1 = sln1
x2, y2 = sln2
plot(x1, y1, label="$u_1$")
plot(x2, y2, label="$u_2$")
legend()
show()
