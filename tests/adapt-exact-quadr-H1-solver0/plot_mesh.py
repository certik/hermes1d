from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("mesh_ref.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="mesh_ref")
data = numpy.loadtxt("mesh.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="mesh")
legend()
show()
