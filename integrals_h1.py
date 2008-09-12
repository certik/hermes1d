def int_grad_u_grad_v(u, v):
    return (u.diff()*v.diff()).integrate()

def int_u_v(u, v):
    return (u*v).integrate()
