from math import sin

from hermes1d.assembly import integ_quad, test3

def f(x):
    return sin(x)

print test3()
print integ_quad(f, 0, 3.1415926)
