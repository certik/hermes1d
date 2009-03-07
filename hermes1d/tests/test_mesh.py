from hermes1d import Node, Element, Mesh

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
    m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2
    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == 3

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
    m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[0].dofs[2] == 4
    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2
    assert m.elements[1].dofs[2] == 5
    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == 3
    assert m.elements[2].dofs[2] == 6

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
    m.assign_dofs()
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
    m.assign_dofs()
    assert m.elements[0].dofs[0] == 0
    assert m.elements[0].dofs[1] == 1
    assert m.elements[0].dofs[2] == 4
    assert m.elements[0].dofs[3] == 5

    assert m.elements[1].dofs[0] == 1
    assert m.elements[1].dofs[1] == 2

    assert m.elements[2].dofs[0] == 2
    assert m.elements[2].dofs[1] == 3
    assert m.elements[2].dofs[2] == 6

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
    m.assign_dofs()
    assert m.elements[0].dofs[0] == -1
    assert m.elements[0].dofs[1] == 0
    assert m.elements[0].dofs[2] == 3
    assert m.elements[0].dofs[3] == 4

    assert m.elements[1].dofs[0] == 0
    assert m.elements[1].dofs[1] == 1

    assert m.elements[2].dofs[0] == 1
    assert m.elements[2].dofs[1] == 2
    assert m.elements[2].dofs[2] == 5
