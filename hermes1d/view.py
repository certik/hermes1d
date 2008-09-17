class View(object):
    """
    Provides basic methods for plotting meshes, functions, basefunctions, ...

    The subclasses of View provide customized behavior, like looping through
    base functions using arrow keys etc.
    """

    def __init__(self, title="View", x=-1, y=-1, width=500, height=500):
        from matplotlib.pyplot import figure
        self.fig = figure()
        self.fig.canvas.manager.resize(width, height)
        self.fig.canvas.manager.set_window_title(title)
        if x != -1 and y != -1:
            # this only works for the GTK backend of matplotlib:
            self.fig.canvas.manager.window.move(x, y)

        def click(canvas, event):
            key = canvas._get_key(event)
            self.key_press_event(key)
        self.fig.canvas.connect("key_press_event", click)

    def key_press_event(self, key):
        if key == "q":
            # comment the following two lines if you just want to close the
            # window:
            import sys
            sys.exit()
            import pylab
            pylab.close(self.fig)

    def plot_mesh(self, mesh):
        x = mesh.get_nodes_x()
        y = [0]*len(x)
        self.fig.gca().plot(x, y, "ko")

    def plot_function(self, f):
        x, y = f.get_xy()
        self.fig.gca().plot(x, y)

    def xplot_base_function(self, space, idx):
        from numpy import zeros
        from hermes1d import Solution
        x = zeros(space.get_num_dofs())
        x[idx] = 1.0
        sln = Solution()
        sln.set_fe_solution(space, x)
        self.plot_function(sln)

class ScalarView(View):

    def __init__(self, title="ScalarView", x=-1, y=-1, width=500, height=500):
        View.__init__(self, title, x, y, width, height)

    def show(self, f):
        self.plot_function(f)
        self.plot_mesh(f.mesh)
        self.fig.canvas.draw()

class BaseView(ScalarView):

    def __init__(self, title="BaseView", x=-1, y=-1, width=500, height=500):
        ScalarView.__init__(self, title, x, y, width, height)

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
        self.fig.clf()
        self.plot_function(self.space.get_base_function(self.base_index))
        self.plot_mesh(self.space.mesh)
        self.fig.gca().set_title("dof = %d" % self.base_index)
        self.fig.gca().set_xlim(self.space.mesh.get_min_max())
        self.fig.gca().set_ylim((-1, 1))
        self.fig.canvas.draw()

    def show(self, space):
        self.space = space
        self.ndofs = space.get_num_dofs()
        self.base_index = 0
        self.update_solution()

class MeshView(View):

    def __init__(self, title="BaseView", x=-1, y=-1, width=500, height=500):
        View.__init__(self, title, x, y, width, height)

    def show(self, mesh):
        self.plot_mesh(mesh)
        self.fig.canvas.draw()

class MatrixView(View):

    def __init__(self, title="MatrixView", x=-1, y=-1, width=500, height=500):
        View.__init__(self, title, x, y, width, height)

    def show(self, dp):
        self.fig.gca().spy(dp.A)
