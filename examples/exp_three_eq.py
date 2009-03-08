from hermes1d import Node, Element, Mesh, DiscreteProblem

from numpy import zeros
a = 0.
b = 1.
N = 40
nodes = [Node((b-a)/N * i) for i in range(N)]
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N-1)]
m1 = Mesh(nodes, elements)
m1.set_bc(left=True, value=1)
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N-1)]
m2 = Mesh(nodes, elements)
m2.set_bc(left=True, value=2)
elements = [Element(nodes[i], nodes[i+1], order=1) for i in range(N-1)]
m3 = Mesh(nodes, elements)
m3.set_bc(left=True, value=-1)

d = DiscreteProblem(meshes=[m1, m2, m3])
k = 1.0
def J(i, j):
    def f11(y1, y2, y3, t):
        return 1
    def f12(y1, y2, y3, t):
        return 0
    def f13(y1, y2, y3, t):
        return 0
    def f21(y1, y2, y3, t):
        return -2
    def f22(y1, y2, y3, t):
        return 2
    def f23(y1, y2, y3, t):
        return 0
    def f31(y1, y2, y3, t):
        return 1
    def f32(y1, y2, y3, t):
        return 1
    def f33(y1, y2, y3, t):
        return 4
    if i == 0 and j == 0:
        return f11
    elif i == 0 and j == 1:
        return f12
    elif i == 0 and j == 2:
        return f13
    elif i == 1 and j == 0:
        return f21
    elif i == 1 and j == 1:
        return f22
    elif i == 1 and j == 2:
        return f23
    elif i == 2 and j == 0:
        return f31
    elif i == 2 and j == 1:
        return f32
    elif i == 2 and j == 2:
        return f33
    raise ValueError("Wrong i, j (i=%d, j=%d)." % (i, j))
def F(i):
    def f1(u1, u2, u3, t):
        return u1
    def f2(u1, u2, u3, t):
        return -2*u1+2*u2
    def f3(u1, u2, u3, t):
        return u1+u2+4*u3
    if i == 0:
        return f1
    elif i == 1:
        return f2
    elif i == 2:
        return f3
    raise ValueError("Wrong i (i=%d)." % (i))
d.set_rhs(F, J)
d.assign_dofs()
J = d.assemble_J()
Y = zeros((J.shape[0],))
error = 1e10
i = 0
while error > 1e-4:
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
sln1, sln2, sln3 = d.linearize(x, 5)
x1, y1 = sln1
x2, y2 = sln2
x3, y3 = sln3
plot(x1, y1, label="$u_1$")
plot(x2, y2, label="$u_2$")
plot(x3, y3, label="$u_3$")
legend()
show()
