from hermes1d import Node, Element, Mesh, DiscreteProblem

def test_node1():
    n = Node(3)
    assert n.x == 3

def test_element1():
    n1 = Node(1)
    n2 = Node(3)
    e = Element(n1, n2)
    assert e.order == 1
    assert e.nodes == (n1, n2)

def test_element2():
    n1 = Node(1)
    n2 = Node(3)
    e = Element(n1, n2, order=2)
    assert e.order == 2
    assert e.nodes == (n1, n2)

def test_mesh1():
    n0 = Node(1)
    n1 = Node(3)
    n2 = Node(4)
    n3 = Node(5)
    e0 = Element(n0, n1)
    e1 = Element(n1, n2)
    e2 = Element(n2, n3)
    nodes = (n0, n1, n2, n3)
    elements = (e0, e1, e2)
    m = Mesh(nodes, elements)
    assert m.nodes == nodes
    assert m.elements == elements
    assert m.nodes[0] == n0
    assert m.nodes[3] == n3

def test_mesh2():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=1)
    e2 = Element(n2, n3, order=1)
    e3 = Element(n3, n4, order=1)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    ndofs = m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2
    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == 3

    assert ndofs == 4

def test_mesh3():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=2)
    e2 = Element(n2, n3, order=2)
    e3 = Element(n3, n4, order=2)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    ndofs = m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[0].dofs[2] == 4
    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2
    assert m.elements[1].dofs[2] == 5
    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == 3
    assert m.elements[2].dofs[2] == 6

    assert ndofs == 7

def test_mesh4():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=3)
    e2 = Element(n2, n3, order=3)
    e3 = Element(n3, n4, order=3)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    ndofs = m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[0].dofs[2] == 4
    assert m.elements[0].dofs[3] == 5

    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2
    assert m.elements[1].dofs[2] == 6
    assert m.elements[1].dofs[3] == 7

    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == 3
    assert m.elements[2].dofs[2] == 8
    assert m.elements[2].dofs[3] == 9

    assert ndofs == 10

def test_mesh5():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=3)
    e2 = Element(n2, n3, order=1)
    e3 = Element(n3, n4, order=2)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    ndofs = m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[0].dofs[2] == 4
    assert m.elements[0].dofs[3] == 5

    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2

    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == 3
    assert m.elements[2].dofs[2] == 6

    assert ndofs == 7

def test_mesh6():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=3)
    e2 = Element(n2, n3, order=1)
    e3 = Element(n3, n4, order=2)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    m.set_bc(left=True, value=1)
    ndofs = m.assign_dofs()
    assert m.elements[0].dofs[0] == -1
    assert m.elements[0].dofs[1] == 0
    assert m.elements[0].dofs[2] == 3
    assert m.elements[0].dofs[3] == 4

    assert m.elements[1].dofs[0] == 0
    assert m.elements[1].dofs[1] == 1

    assert m.elements[2].dofs[0] == 1
    assert m.elements[2].dofs[1] == 2
    assert m.elements[2].dofs[2] == 5

    assert ndofs == 6

def test_mesh7():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=3)
    e2 = Element(n2, n3, order=1)
    e3 = Element(n3, n4, order=2)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    m.set_bc(left=False, value=1)
    ndofs = m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[0].dofs[2] == 3
    assert m.elements[0].dofs[3] == 4

    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2

    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == -1
    assert m.elements[2].dofs[2] == 5

    assert ndofs == 6

def test_mesh8():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=3)
    e2 = Element(n2, n3, order=1)
    e3 = Element(n3, n4, order=2)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m1 = Mesh(nodes, elements)
    m1.set_bc(left=False, value=1)
    e4 = Element(n1, n2, order=3)
    e5 = Element(n2, n3, order=1)
    e6 = Element(n3, n4, order=2)
    elements = (e4, e5, e6)
    m2 = Mesh(nodes, elements)
    m2.set_bc(left=True, value=1)

    m1.assign_dofs()
    m2.assign_dofs()

    assert m1.elements[0].dofs[0] == 0
    assert m1.elements[0].dofs[1] == 1
    assert m1.elements[0].dofs[2] == 3
    assert m1.elements[0].dofs[3] == 4

    assert m1.elements[1].dofs[0] == 1
    assert m1.elements[1].dofs[1] == 2

    assert m1.elements[2].dofs[0] == 2
    assert m1.elements[2].dofs[1] == -1
    assert m1.elements[2].dofs[2] == 5

    assert m2.elements[0].dofs[0] == -1
    assert m2.elements[0].dofs[1] == 0
    assert m2.elements[0].dofs[2] == 3
    assert m2.elements[0].dofs[3] == 4

    assert m2.elements[1].dofs[0] == 0
    assert m2.elements[1].dofs[1] == 1

    assert m2.elements[2].dofs[0] == 1
    assert m2.elements[2].dofs[1] == 2
    assert m2.elements[2].dofs[2] == 5

def test_mesh9():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=3)
    e2 = Element(n2, n3, order=1)
    e3 = Element(n3, n4, order=2)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m1 = Mesh(nodes, elements)
    m1.set_bc(left=False, value=1)
    e4 = Element(n1, n2, order=3)
    e5 = Element(n2, n3, order=1)
    e6 = Element(n3, n4, order=2)
    elements = (e4, e5, e6)
    m2 = Mesh(nodes, elements)
    m2.set_bc(left=True, value=1)

    d = DiscreteProblem(meshes=[m1, m2])
    ndofs = d.assign_dofs()

    assert m1.elements[0].dofs[0] == 0
    assert m1.elements[0].dofs[1] == 1
    assert m1.elements[0].dofs[2] == 3
    assert m1.elements[0].dofs[3] == 4

    assert m1.elements[1].dofs[0] == 1
    assert m1.elements[1].dofs[1] == 2

    assert m1.elements[2].dofs[0] == 2
    assert m1.elements[2].dofs[1] == -1
    assert m1.elements[2].dofs[2] == 5

    assert m2.elements[0].dofs[0] == -1
    assert m2.elements[0].dofs[1] == 0 + 6
    assert m2.elements[0].dofs[2] == 3 + 6
    assert m2.elements[0].dofs[3] == 4 + 6

    assert m2.elements[1].dofs[0] == 0 + 6
    assert m2.elements[1].dofs[1] == 1 + 6

    assert m2.elements[2].dofs[0] == 1 + 6
    assert m2.elements[2].dofs[1] == 2 + 6
    assert m2.elements[2].dofs[2] == 5 + 6

    assert ndofs == 12
    assert d.get_mesh_number(0) == 0
    assert d.get_mesh_number(4) == 0
    assert d.get_mesh_number(5) == 0
    assert d.get_mesh_number(6) == 1
    assert d.get_mesh_number(11) == 1

def test_discrete_problem():
    n1 = Node(0)
    n2 = Node(1)
    n3 = Node(2)
    n4 = Node(3)
    e1 = Element(n1, n2, order=1)
    e2 = Element(n2, n3, order=1)
    e3 = Element(n3, n4, order=1)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m1 = Mesh(nodes, elements)
    m1.set_bc(left=True, value=0)
    e4 = Element(n1, n2, order=1)
    e5 = Element(n2, n3, order=1)
    e6 = Element(n3, n4, order=1)
    elements = (e4, e5, e6)
    m2 = Mesh(nodes, elements)
    m2.set_bc(left=True, value=1)

    d = DiscreteProblem(meshes=[m1, m2])
    def J(i, j):
        def f11(y1, y2, t):
            return 0
        def f12(y1, y2, t):
            return 1
        def f21(y1, y2, t):
            return -k**2
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
    def F(i):
        def f1(y1, y2, t):
            return y2
        def f2(y1, y2, t):
            k = 2.0
            return -k**2 * y1
        if i == 0:
            return f1
        elif i == 1:
            return f2
        raise ValueError("Wrong i (i=%d)." % (i))
    d.set_rhs(F, J)
    d.assign_dofs()
    J = d.assemble_J()
    F = d.assemble_F()
    # solve: J*dY = -F
    print J
    print F
