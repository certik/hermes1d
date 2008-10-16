def convert_mat(mtx):
    """
    Converts a scipy matrix "mtx" to a pysparse matrix.
    """
    from pysparse import spmatrix
    mtx = mtx.tocsr()
    A = spmatrix.ll_mat(*mtx.shape)
    for i in xrange( mtx.indptr.shape[0] - 1 ):
        ii = slice( mtx.indptr[i], mtx.indptr[i+1] )
        n_in_row = ii.stop - ii.start
        A.update_add_at( mtx.data[ii], [i] * n_in_row, mtx.indices[ii] )
    return A

def eigen(A, B):
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
    tau = -40
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
    E_exact = [-1./2/(n-0.5)**2/4 for n in [1]+[2]*3+[3]*5 + [4]*8 + [5]*15]
    print "eigenvalues (i, FEM, exact, error):"
    for i, E in enumerate(lmbd):
        a = E
        b = E_exact[i]
        print "%2d: %10f %10f %f%%" % (i, a, b, abs((a-b)/b)*100)
    return Q
