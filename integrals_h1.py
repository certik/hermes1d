def int_grad_u_grad_v(u, v):
    return (u.diff()*v.diff()).integrate()

def int_u_v(u, v):
    return (u*v).integrate()

def int_F_u_v(F, u, v):
    return (F*u*v).integrate()

def int_grad_u_v_overx(u, v):
    from hermes1d import LinearFunction
    x = LinearFunction(u.mesh)
    return (u*v/x).integrate()
