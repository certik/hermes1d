class ScalarView:

    def __init__(self, name):
        pass

    def show(self, f):
        from numpy import arange
        import pylab
        x, y = f.get_xy()
        pylab.plot(x, y)
        def click(event):
            if event.key == "q":
                import sys
                sys.exit()
        pylab.connect("key_press_event", click)
        pylab.show()
