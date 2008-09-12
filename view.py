class ScalarView:

    def __init__(self, name):
        pass

    def show(self, f):
        from numpy import arange
        import pylab
        a, b = f.domain_range()
        x = arange(a, b, float(b-a)/100)
        y = [f(_x) for _x in x]
        pylab.plot(x, y)
        def click(event):
            if event.key == "q":
                import sys
                sys.exit()
        pylab.connect("key_press_event", click)
        pylab.show()
