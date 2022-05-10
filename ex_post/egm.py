# Import package and module
import numpy as np
import utility as util
import tools


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
