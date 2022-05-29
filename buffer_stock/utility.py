from numba import njit

@njit(fastmath=True)
def func(c,par):
    return c**(1-par.rho)/(1-par.rho)

@njit(fastmath=True)
def marg_util(c,par):
    return c**(-par.rho)

@njit(fastmath=True)
def inv_marg_util(u,par):
    return u**(-1/par.rho)
