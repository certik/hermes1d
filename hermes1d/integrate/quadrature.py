## Automatically adapted for scipy Oct 21, 2005 by

# Author: Travis Oliphant

from scipy.special.orthogonal import p_roots as p_roots_orig
from numpy import sum, isinf, isscalar, asarray, real, empty

_cache = {}
@profile
def p_roots(n):
    try:
        return _cache[n]
    except KeyError:
        _cache[n] = p_roots_orig(n)
        return _cache[n]

@profile
def fixed_quad(func,a,b,args=(),n=5):
    """Compute a definite integral using fixed-order Gaussian quadrature.

  Description:

    Integrate func from a to b using Gaussian quadrature of order n.

  Inputs:

    func -- a Python function or method to integrate
            (must accept vector inputs)
    a -- lower limit of integration
    b -- upper limit of integration
    args -- extra arguments to pass to function.
    n -- order of quadrature integration.

  Outputs: (val, None)

    val -- Gaussian quadrature approximation to the integral.

  See also:
    
    quad - adaptive quadrature using QUADPACK
    dblquad, tplquad - double and triple integrals
    romberg - adaptive Romberg quadrature
    quadrature - adaptive Gaussian quadrature
    romb, simps, trapz - integrators for sampled data
    cumtrapz - cumulative integration for sampled data
    ode, odeint - ODE integrators
    """
    [x,w] = p_roots(n)
    x = real(x)
    ainf, binf = map(isinf,(a,b))
    if ainf or binf:
        raise ValueError, "Gaussian quadrature is only available for " \
              "finite limits."
    y = (b-a)*(x+1)/2.0 + a
    return (b-a)/2.0*sum(w*func(y,*args),0), None

def vectorize1(func, args=(), vec_func=False):
    if vec_func:
        def vfunc(x):
            return func(x, *args)
    else:
        def vfunc(x):
            if isscalar(x):
                return func(x, *args)
            x = asarray(x)
            # call with first point to get output type
            y0 = func(x[0], *args)
            n = len(x)
            output = empty((n,), dtype=y0.dtype)
            output[0] = y0
            for i in xrange(1, n):
                output[i] = func(x[i], *args)
            return output
    return vfunc

@profile
def quadrature(func,a,b,args=(),tol=1.49e-8,maxiter=50, vec_func=True):
    """Compute a definite integral using fixed-tolerance Gaussian quadrature.

  Description:

    Integrate func from a to b using Gaussian quadrature
    with absolute tolerance tol.

  Inputs:

    func -- a Python function or method to integrate.
    a -- lower limit of integration.
    b -- upper limit of integration.
    args -- extra arguments to pass to function.
    tol -- iteration stops when error between last two iterates is less than
           tolerance.
    maxiter -- maximum number of iterations.
    vec_func -- True or False if func handles arrays as arguments (is
                a "vector" function ). Default is True.

  Outputs: (val, err)

    val -- Gaussian quadrature approximation (within tolerance) to integral.
    err -- Difference between last two estimates of the integral.

  See also:
    
    romberg - adaptive Romberg quadrature
    fixed_quad - fixed-order Gaussian quadrature
    quad - adaptive quadrature using QUADPACK
    dblquad, tplquad - double and triple integrals
    romb, simps, trapz - integrators for sampled data
    cumtrapz - cumulative integration for sampled data
    ode, odeint - ODE integrators
    """
    err = 100.0
    val = err
    n = 1
    vfunc = vectorize1(func, args, vec_func=vec_func)
    while (err > tol) and (n < maxiter):
        newval = fixed_quad(vfunc, a, b, (), n)[0]
        err = abs(newval-val)
        val = newval
        n = n + 1
    if n == maxiter:
        print "maxiter (%d) exceeded. Latest difference = %e" % (n,err)
    return val, err
