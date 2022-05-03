from numba import njit

@njit
def func(c,par):
    return c**(1-par.rho)/(1-par.rho)

def marg_util(c,par):
    return c**(-par.rho)

def inv_marg_util(u,par):
    return u**(-1/par.rho)
