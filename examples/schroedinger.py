"""
Solves the second order ODE:

y'' + y = 0
y(0)=0; y'(0)=1

So the solution is y(x) = sin(x)
"""
from hermes1d import Node, Element, Mesh, DiscreteProblem

from math import pi
from numpy import zeros
from numpy.linalg import inv, eig, eigh

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

# enumeration of unknowns:
d.assign_dofs()

print "assembling"
A = d.assemble_schroed(rhs=False)
B = d.assemble_schroed(rhs=True)
print "inverting"
M = inv(B)*A
print "solving:"
w, v = eigh(M)
print w
print v
