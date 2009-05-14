from numpy import matrix, array, dot
from numpy.linalg import inv, eig

from hermes1d import Node, Element, DiscreteProblem, Mesh

def eq(a, b, eps=1e-8, verbose=False):
    if verbose:
        print a, b, abs(a-b)
    return abs(a-b) < eps

def test_hydrogen1():
    a = 0.
    b = 30
    N = 20
    x_values =[(b-a)/N * i for i in range(N+1)]
    nodes = [Node(x) for x in x_values]
    elements = [Element(nodes[i], nodes[i+1], order=2) for i in range(N)]
    m1 = Mesh(nodes, elements)

    def schroed_l(m, l=0):
        d = DiscreteProblem(meshes=[m])
        d.assign_dofs()
        A = d.assemble_schroed(rhs=False, l=l, pot="hydrogen", a=b)
        B = d.assemble_schroed(rhs=True, l=l, pot="hydrogen", a=b)
        M = dot(inv(B), matrix(A))
        w, v = eig(M)
        r = list(w)
        return r

    r = schroed_l(m1, l=0)
    r.extend(schroed_l(m1, l=1))
    r.extend(schroed_l(m1, l=2))
    r.extend(schroed_l(m1, l=3))
    r.sort(key=lambda x: x.real)
    E = [-1./(2*n**2) for n in [1] + [2]*2 + [3]*3 + [4]*4]
    for i in range(len(E)):
        assert eq(r[i], E[i], eps=0.01)

def test_hydrogen2():
    a = 0.
    b = 30
    N = 20
    x_values =[(b-a)/N * i for i in range(N+1)]
    nodes = [Node(x) for x in x_values]
    elements = [Element(nodes[i], nodes[i+1], order=2) for i in range(N)]
    m1 = Mesh(nodes, elements)

    def schroed_l(m, l=0):
        m.set_bc(left=False, value=0)
        d = DiscreteProblem(meshes=[m])
        d.assign_dofs()
        A = d.assemble_schroed(rhs=False, l=l, pot="hydrogen", a=b,
                bc_calculate=False)
        B = d.assemble_schroed(rhs=True, l=l, pot="hydrogen", a=b,
                bc_calculate=False)
        M = dot(inv(B), matrix(A))
        w, v = eig(M)
        r = list(w)
        return r

    r = schroed_l(m1, l=0)
    r.extend(schroed_l(m1, l=1))
    r.extend(schroed_l(m1, l=2))
    r.extend(schroed_l(m1, l=3))
    r.sort(key=lambda x: x.real)
    E = [-1./(2*n**2) for n in [1] + [2]*2 + [3]*3 + [4]*4]
    for i in range(len(E)):
        assert eq(r[i], E[i], eps=0.01, verbose=True)

def test_hydrogen3():
    a = 0.
    b = 30
    N = 20
    x_values =[(b-a)/N * i for i in range(N+1)]
    nodes = [Node(x) for x in x_values]
    elements = [Element(nodes[i], nodes[i+1], order=2) for i in range(N)]
    m1 = Mesh(nodes, elements)

    def schroed_l(m, l=0):
        if l == 0:
            m.set_bc(left=True, value=-1)
        else:
            m.set_bc(left=True, value=0)
        m.set_bc(left=False, value=0)
        d = DiscreteProblem(meshes=[m])
        d.assign_dofs()
        A = d.assemble_schroed(rhs=False, l=l, pot="hydrogen", a=b,
                bc_calculate=False)
        B = d.assemble_schroed(rhs=True, l=l, pot="hydrogen", a=b,
                bc_calculate=False)
        M = dot(inv(B), matrix(A))
        w, v = eig(M)
        r = list(w)
        return r

    r = schroed_l(m1, l=0)
    r.extend(schroed_l(m1, l=1))
    r.extend(schroed_l(m1, l=2))
    r.extend(schroed_l(m1, l=3))
    r.sort(key=lambda x: x.real)
    E = [-1./(2*n**2) for n in [1] + [2]*2 + [3]*3 + [4]*4]
    for i in range(len(E)):
        assert eq(r[i], E[i], eps=0.17, verbose=True)
