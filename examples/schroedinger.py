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
b = pi
b = 5.

# number of elements:
N = 5

# x values of the nodes:
x_values =[(b-a)/N * i for i in range(N+1)]

# define nodes:
nodes = [Node(x) for x in x_values]

# define elements of the 1st mesh
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N)]
m1 = Mesh(nodes, elements)

def schroed_l(m, l=0):
    if l == 0:
        m.set_bc(left=True, value=0)
    else:
        m.set_bc(left=True, value=0)
    m.set_bc(left=False, value=0)

    # definition of the ODE system:
    d = DiscreteProblem(meshes=[m])

    # enumeration of unknowns:
    d.assign_dofs()
    #print m
    m.print_dofs()
    #stop

    print "assembling"
    A = d.assemble_schroed(rhs=False, pot="well1d")
    B = d.assemble_schroed(rhs=True, pot="well1d")
    print A
    print B
    #stop
    print "inverting"
    from numpy import matrix, array
    M = matrix(inv(B))*matrix(A)
    M = array(M)
    print "solving:"
    w, v = eigh(M)
    print w
    vec = v[:, 0]
    print vec
    print "test"
    #import pdb
    #pdb.set_trace()
    print matrix(A)*matrix(vec).T
    print w[0]*matrix(B)*matrix(vec).T
    #stop
    # sort w and v:
    r = []
    for i in range(len(w)):
        #vec = v[i, :]
        vec = v[:, i]
        #print "ooo"
        #print v
        #print v[:, i]
        #stop
        x, y = d.linearize(vec, 10)[0]
        r.append((w[i], vec, l, x, y))
    r.sort(key=lambda x: x[0])
    return r

r = schroed_l(m1, l=0)
#r.extend(schroed_l(m1, l=1))
#r.extend(schroed_l(m1, l=2))
#r.extend(schroed_l(m1, l=3))
r.sort(key=lambda x: x[0])
print "results:"
for i in range(4):
    w, v, l, x, y = r[i]
    print "l=%d; E=%f" % (l, w)
print "plotting:"
from pylab import plot, show, legend
for i in range(1):
    w, v, l, x, y = r[i]
    plot(x, y, label="l=%d, eig=%f" % (l, w))
legend()
show()
