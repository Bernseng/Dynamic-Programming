import numpy as np
from numba import njit

@njit(fastmath=True)
def p_plus_func(p,psi,par,t):
    if t+1<=par.Tr:
        p_plus = p*psi*par.G*par.L[t]
        p_plus = np.fmax(p_plus,par.p_min) # lower bound
        p_plus = np.fmin(p_plus,par.p_max) # upper bound
    else: 
        p_plus = p*par.G*par.L[t] # no shocks to permanent income after retirement
        p_plus = np.fmax(p_plus,par.p_min) # lower bound
        p_plus = np.fmin(p_plus,par.p_max) # upper bound
    return p_plus

@njit(fastmath=True)
def n_plus_func(d,par,z):
    n_plus = (par.Rh)*(1-par.delta)*d   # deleted "+z"
    n_plus = np.fmin(n_plus,par.n_max) # upper bound
    return n_plus

@njit(fastmath=True)
def m_plus_func(a,p_plus,xi_plus,par):
    y_plus = p_plus*xi_plus
    m_plus = par.R*a+ y_plus
    return m_plus

@njit(fastmath=True)
def x_plus_func(m_plus,n_plus,par):
    return m_plus + (1-par.tau)*n_plus