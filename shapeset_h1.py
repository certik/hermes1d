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
            elif idx == 5:
                return (x**2-1)*(7*x**2-3)*sqrt(9./2)/8
            elif idx == 6:
                return (x**2-1)*(21*x**4-14*x**2+1)*sqrt(11./2)/16
            elif idx == 7:
                return (x**2-1)*(33*x**4-30*x**2+5)*sqrt(13./2)/16
            elif idx == 8:
                return (x**2-1)*(429*x**6-495*x**4+135*x**2-5)*sqrt(15./2)/128
            elif idx == 9:
                return (x**2-1)*(715*x**6-1001*x**4+385*x**2-35)*sqrt(17./2)/128
            elif idx == 10:
                return (x**2-1)*(2431*x**8-4004*x**6+2002*x**4-308*x**2+7)*sqrt(19./2)/256
            else:
                raise NotImplementedError("Such idx is not implemented (yet).")
        elif diff == 1:
            if idx == 0:
                return -0.5
            elif idx == 1:
                return 0.5
            else:
                raise NotImplementedError("Such idx is not implemented (yet).")
        elif diff > 1:
            return 0.
        raise NotImplementedError()
