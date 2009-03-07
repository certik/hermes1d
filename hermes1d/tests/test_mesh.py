from hermes1d import Node, Element

def test_node1():
    n = Node(3)
    assert n.x == 3

def test_element1():
    n1 = Node(1)
    n2 = Node(3)
    e = Element(n1, n2)
    assert e.nodes == (n1, n2)
