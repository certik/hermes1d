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
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2)
    e2 = Element(n2, n3)
    e3 = Element(n3, n4)
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    assert m.nodes == nodes
    assert m.elements == elements

def test_mesh2():
    n1 = Node(1)
    n2 = Node(3)
    n3 = Node(4)
    n4 = Node(5)
    e1 = Element(n1, n2, order=2)
    e2 = Element(n2, n3, order=2)
    e3 = Element(n3, n4, order=2)
    assert e3.order == 2
    nodes = (n1, n2, n3, n4)
    elements = (e1, e2, e3)
    m = Mesh(nodes, elements)
    assert m.nodes == nodes
    assert m.elements == elements
