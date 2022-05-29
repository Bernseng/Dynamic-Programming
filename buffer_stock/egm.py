# Import package and module
import numpy as np
from numba import njit, prange

import utility 
import trans
import tools

from consav import linear_interp

"""
def EGM(sol,t,par):
    if t+1 <= par.Tr: 
        fac = np.tile(par.G*par.L[t]*par.psi_vec, par.Na) 
        xi = np.tile(par.xi_vec,par.Na)
        a = np.repeat(par.grid_a[t],par.Nshocks) 

        w = np.tile(par.w,(par.Na,1))
        dim = par.Nshocks
    else:
        fac = par.G*par.L[t]*np.ones((par.Na))
        xi = np.ones((par.Na))
        a = par.grid_a[t,:]
            
        w = np.ones((par.Na,1))
        dim = 1

    inv_fac = 1/fac

    # Futute m and c
    m_plus = inv_fac*par.R*a+xi
    c_plus = tools.interp_linear_1d(sol.m[t+1,:],sol.c[t+1,:], m_plus)

    # Future marginal utility
    marg_u_plus = util.marg_util(fac*c_plus,par)
    marg_u_plus = np.reshape(marg_u_plus,(par.Na,dim))
    avg_marg_u_plus = np.sum(w*marg_u_plus,1)

    # Currect C and m
    sol.c[t,1:]= util.inv_marg_util(par.beta*par.R*avg_marg_u_plus,par)
    sol.m[t,1:]=par.grid_a[t,:]+sol.c[t,1:]

    # f. current v
    # if sol.c[t,1:] > 0:
    #     sol.inv_v[t,1:] = 1.0/(util.func(c[i_a],par) + par.beta*avg_v_plus)
    # else:
    #     sol.inv_v[t,1:] = 0

    return sol
"""

@njit(parallel=True)
def egm(par,sol,t,m,c,inv_v):
    """ apply egm step """
        
    # loop over end-of-period assets
    for i_a in prange(par.Na): # parallel

        a = par.grid_a[t,i_a]
        still_working_next_period = t+1 <= par.Tr-1
        Nshocks = par.Nshocks if still_working_next_period else 1

        # loop over shocks
        avg_marg_u_plus = 0
        avg_v_plus = 0
        for i_shock in range(Nshocks):
            
            # a. prep
            if still_working_next_period:
                fac = par.G*par.L[t]*par.psi[i_shock]
                w = par.w[i_shock]
            else:
                fac = par.G*par.L[t]
                w = 1
        
            psi_plus = par.psi[i_shock]
            xi_plus = par.xi[i_shock]

            # b. future m and c
            m_plus = trans.m_plus_func(a,xi_plus,psi_plus,par,t)
            c_plus = linear_interp.interp_1d(sol.m[t+1,:],sol.c[t+1,:],m_plus)
            inv_v_plus = linear_interp.interp_1d(sol.m[t+1,:],sol.inv_v[t+1,:],m_plus)
            v_plus = 1.0/inv_v_plus
            
            # c. average future marginal utility
            marg_u_plus = utility.marg_util(fac*c_plus,par)
            avg_marg_u_plus += w*marg_u_plus
            avg_v_plus += w*(fac**(1-par.rho))*v_plus

        # d. current c
        c[i_a] = utility.inv_marg_util(par.beta*par.R*avg_marg_u_plus,par)

        # e. current m
        m[i_a] = a + c[i_a]

        # f. current v
        if c[i_a] > 0:
            inv_v[i_a] = 1.0/(utility.func(c[i_a],par) + par.beta*avg_v_plus)
        else:
            inv_v[i_a] = 0
