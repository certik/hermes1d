class View(object):
    pass

class ScalarView(View):

    def __init__(self, *args):
        from matplotlib.pyplot import figure
        self.fig = figure()
        def click(canvas, event):
            key = canvas._get_key(event)
            self.key_press_event(key)
        self.fig.canvas.connect("key_press_event", click)

    def key_press_event(self, key):
        if key == "q":
            import pylab
            pylab.close(self.fig)

    def plot(self, f):
        x, y = f.get_xy()
        self.fig.gca().plot(x, y)

    def interact(self):
        pass


    def show(self, f):
        self.plot(f)
        self.interact()

class BaseView(ScalarView):

    def __init__(self, *args):
        ScalarView.__init__(self, *args)

    def key_press_event(self, key):
        if key == "left":
            if self.base_index > 0:
                self.base_index -=1
                self.update_solution()
        elif key == "right":
            if self.base_index < self.ndofs - 1:
                self.base_index +=1
                self.update_solution()
        else:
            ScalarView.key_press_event(self, key)

    def update_solution(self):
        from numpy import zeros
        x = zeros(self.ndofs)
        x[self.base_index] = 1.0
        self.sln.set_fe_solution(self.space, x)
        self.fig.clf()
        self.plot(self.sln)
        self.fig.gca().set_title("dof = %d" % self.base_index)
        self.fig.canvas.draw()

    def show(self, space):
        from hermes1d import Solution
        self.space = space
        self.ndofs = space.get_num_dofs()
        self.base_index = 0
        self.sln = Solution()
        self.update_solution()
        ScalarView.interact(self)
