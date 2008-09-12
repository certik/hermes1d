class Function(object):
    pass

class Solution(Function):

    def domain_range(self):
        return 0, 3.14159

    def __call__(self, x):
        idx = int(x)
        if idx < 0:
            idx = 0
        if idx > len(self.coeff) - 1:
            idx = len(self.coeff) - 1
        return self.coeff[idx]

    def set_fe_solution(self, space, x):
        self.coeff = x

class ShapeFunction(Function):

    def __init__(self, shapeset, idx):
        self.idx = idx
