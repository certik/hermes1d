class Solution(object):

    def domain_range(self):
        return 0, 3.14159

    def __call__(self, x):
        from math import sin
        return sin(x)
