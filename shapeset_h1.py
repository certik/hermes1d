from math import sqrt

class H1Shapeset(object):
    """
    Represents shape functions defined on the reference element.
    """

    def get_value_reference(self, x, idx, diff=0):
        """
        Returns the value (diff=0) of a shapefunction on a reference element.

        x is a coordinate on a reference element, i.e. -1 <= x <= 1.

        For diff>0, returns the value of derivatives.
        """
        assert -1 <= x <= 1
        if diff == 0:
            if idx == 0:
                return (1-x)/2
            elif idx == 1:
                return (1+x)/2
            elif idx == 2:
                return (x**2-1)*sqrt(3./2)/2
            elif idx == 3:
                return (x**2-1)*x*sqrt(5./2)/2
            elif idx == 4:
                return (x**2-1)*(5*x**2-1)*sqrt(7./2)/8
            else:
                raise NotImplementedError("Such idx is not implemented (yet).")
        elif diff == 1:
            if idx == 0:
                return -1
            if idx == 1:
                return 1
        elif diff > 1:
            return 0.
        raise NotImplementedError()
