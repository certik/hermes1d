#! /usr/bin/env python

from hermes1d import finalize, Mesh, H1Shapeset, H1Space, DiscreteProblem, \
        Solution, ScalarView, int_grad_u_grad_v, int_u_v, BaseView, MeshView, \
        MatrixView, BC_DIRICHLET, BC_NEUMANN, LinearFunction, int_F_u_v, \
        CustomFunction, int_grad_u_v_overx

def dense2sparse(mtx):
    from scipy.sparse import lil_matrix
    m = lil_matrix(shape=mtx.shape, dtype="double")
    n = mtx.shape[0]
    print n
    for i in range(n):
        for j in range(n):
            if abs(mtx[i,j]) > 1e-8:
                m[i,j] = mtx[i,j]
    return m

def convert_mat(mtx):
    """
    Converts a scipy matrix "mtx" to a pysparse matrix.
    """
    from pysparse import spmatrix
    mtx = dense2sparse(mtx)
    mtx = mtx.tocsr()
    A = spmatrix.ll_mat(*mtx.shape)
    for i in xrange( mtx.indptr.shape[0] - 1 ):
        ii = slice( mtx.indptr[i], mtx.indptr[i+1] )
        n_in_row = ii.stop - ii.start
        A.update_add_at( mtx.data[ii], [i] * n_in_row, mtx.indices[ii] )
    return A

def solve(A, B):
    """
    Solves the generalized eigenvalue problem.

    A, B ... scipy matrices
    """
    from pysparse import jdsym, precon, itsolvers
    print "converting to pysparse"
    n = A.shape[0]
    A = convert_mat(A)
    B = convert_mat(B)
    print "solving (%d x %d)" % (n, n)
    Atau = A.copy()
    tau = -1
    Atau.shift(-tau, B)
    K = precon.jacobi(Atau)
    A = A.to_sss()
    B = B.to_sss()
    n_eigs = 4
    kconv, lmbd, Q, it, it_in = jdsym.jdsym(A, B, K, n_eigs, tau, 1e-6, 150,
            itsolvers.qmrs)
    print "number of converged eigenvalues:", kconv
    #levels = []
    #for n1 in range(1, 10):
    #    for n2 in range(1, 10):
    #        levels.append(n1**2 + n2**2)
    #levels.sort()

    # well
    #E_exact = [pi**2/2 * m for m in levels]

    # oscillator
    #E_exact = [1] + [2]*2 + [3]*3 + [4]*4 + [5]*5 + [6]*6

    # hydrogen
    #E_exact = [-1./2/(n-0.5)**2/4 for n in [1]+[2]*3+[3]*5 + [4]*8 + [5]*15]
    E_exact = [-1./2/n**2 for n in [1]+[2]*2**2+[3]*3**3+[4]*4**4+[5]*5**5]
    print "eigenvalues (i, FEM, exact, error):"
    for i, E in enumerate(lmbd):
        a = E
        b = E_exact[i]
        print "%2d: %10f %10f %f%%" % (i, a, b, abs((a-b)/b)*100)
    return Q

def F(r):
    if r == 0:
        return -1.
    V = -1./r
    l = 1.
    return V + l*(l+1)/(2*r**2)

def bilinear_formA(u, v):
    return int_grad_u_grad_v(u, v)/2-int_grad_u_v_overx(u,v)+int_F_u_v(FF, u, v)

def bilinear_formB(u, v):
    return int_u_v(u, v)

print "Creating mesh"
mesh = Mesh()
mesh.set_interval(0, 100)
#mesh.refine_element(0)
mesh.refine_all_elements()
mesh.refine_all_elements()
mesh.refine_all_elements()
#mesh.refine_all_elements()
#mesh.refine_all_elements()
#mesh.refine_all_elements()
shapeset = H1Shapeset()

print "Creating base functions"
space = H1Space(mesh, shapeset)
space.set_uniform_order(1)
space.build()
#bview = BaseView("Base", 0, 0, 400, 400)
#bview.show(space)
#finalize()
FF = CustomFunction(F, space)
#v = ScalarView()
#v.show(FF)
#finalize()


print "Setting up forms"
dpA = DiscreteProblem()
dpA.set_num_equations(1)
dpA.set_spaces(space)
dpA.set_bilinear_form(0, 0, bilinear_formA)
dpA.create_matrix()
print "Assembling"
dpA.assemble_matrix_and_rhs()
print "Setting up forms"
dpB = DiscreteProblem()
dpB.set_num_equations(1)
dpB.set_spaces(space)
dpB.set_bilinear_form(0, 0, bilinear_formB)
dpB.create_matrix()
print "Assembling"
dpB.assemble_matrix_and_rhs()
A = dpA.A
B = dpB.A

print "Solving"
print A
print B
sols = solve(A, B)

s = []
n = sols.shape[1]
for i in range(n):
    sln = Solution()
    vec = sols[:, i]
    sln.set_fe_solution(space, vec)
    s.append(sln)

print "Plotting"
# visualize the solution
w = 400
h = 300
mview = MeshView("Mesh", 0, 0, w, h)
mview.show(mesh)
bview = BaseView("Base", w, 0, w, h)
bview.show(space)
for i, sln in enumerate(s):
    view = ScalarView("Eigenvector %d" % i, i*w, 2*h, w, h)
    view.show(sln)
#view = ScalarView("Error plot", 2*w, h, w, h)
#x = LinearFunction(mesh)
#sln_exact = (1-x**2)/2
#error = (sln-sln_exact)**2
#view.show(error)
#print error.is_zero(1e-9)
#mat1 = MatrixView("Matrix A (%dx%d)" % dpA.A.shape, 0, h, w, h)
#mat1.show(dpA)

finalize()
