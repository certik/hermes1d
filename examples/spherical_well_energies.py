"""
This example calculates the energies for a spherical potential well.

It depends on sympy.
"""
from math import pi

from numpy import array

from sympy import jn_zeros

def E(zeros):
    zeros = array(zeros)
    a = pi  # "a" is the radius of the sphere
    E = zeros**2 / (2*a**2)
    return E

energies = []
for l in range(20):
    for e in E(jn_zeros(l, 4, method="scipy")):
        energies.append((l, e))

energies.sort(key=lambda x:x[1])
for l, e in energies[:13]:
    print "l=%d; E=%f" % (l, e)
