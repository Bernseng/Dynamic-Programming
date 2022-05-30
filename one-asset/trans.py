import numpy as np
from numba import njit

@njit(fastmath=True)
def p_plus_func(p,psi,par,t):
    if t<=par.Tr:
        p_plus = p + np.log(psi) + np.log(par.G) + np.log(par.L[t])
        p_plus = np.fmax(p_plus,par.p_min) #lower bound
        p_plus = np.fmin(p_plus,par.p_max) #upper bound
    else:
        p_plus = p + np.log(par.G) + np.log(par.L[t]) #no shocks to permanent income
        p_plus = np.fmax(p_plus,par.p_min) #lower bound
        p_plus = np.fmin(p_plus,par.p_max) #upper bound
    return p_plus 

@njit(fastmath=True)
def m_plus_func(a,xi_plus,psi_plus,par,t):
    if t<=par.Tr:
        m_plus = par.R*a/(psi_plus*par.G*par.L[t]) + xi_plus
    else:
        m_plus = par.R*a/(par.G*par.L[t]) + 1
    return m_plus

@njit(fastmath=True)
def y_plus_func(p_plus,xi_plus,par,t):
    if t<=par.Tr:    
        y_plus = p_plus + np.log(xi_plus)
    else:
        y_plus = p_plus
    return y_plus

'''
if t < par.simT-1:

                if t+1 > par.TR-1:
                    m[i,t+1] = par.R*a[i,t] / (par.G*par.L[t]) +  1
                    p[i,t+1] = np.log(par.G) + np.log(par.L[t]) + p[i,t]
                    y[i,t+1] = p[i,t+1]
                else:
                    m[i,t+1] = par.R*a[i,t] / (par.G*par.L[t]*sim.psi[i,t+1]) + sim.xi[i,t+1]
                    p[i,t+1] = np.log(par.G) + np.log(par.L[t]) + p[i,t] + np.log(sim.psi[i,t+1])   
                    if sim.xi[i,t+1] > 0:
                        y[i,t+1] = p[i,t+1] + np.log(sim.xi[i,t+1])
                    else:
                        y[i,t+1] = -np.inf
'''