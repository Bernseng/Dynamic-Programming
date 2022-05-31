import numpy as np
from numba import njit,prange

from consav import linear_interp
import tools
import trans

#@njit(parallel=True)
def life_cycle(par,sol,sim):
    
    # unpack 
    m = sim.m
    c = sim.c
    a = sim.a
    p = sim.p
    y = sim.y
    mpc = sim.mpc
    
    sol_c = sol.c
    sol_m = sol.m 
    
    m[0,:] = sim.m0
    p[0,:] = sim.p0
    y[0,:] = sim.p0 + np.log(sim.xi[0,:])
    
    for t in range(par.T):
        
        grid_m = sol_m[t,:]      
        grid_c = sol_c[t,:]

        c[t,:] = tools.interp_linear_1d(grid_m,grid_c,m[t,:])
        a[t,:] = m[t,:] - c[t,:]

        for i in range(par.simN):
            if t<par.T-1:
                m[t+1,i] = trans.m_plus_func(a[t,i],sim.xi[t+1,i],sim.psi[t+1,i],par,t)
                p[t+1,i] = trans.p_plus_func(p[t,i],sim.psi[t+1,i],par,t)
                y[t+1,i] = trans.y_plus_func(p[t+1,i],sim.xi[t+1,i],par,t)
        
        mpc[t,:] = (tools.interp_linear_1d(grid_m,grid_c,m[t,:]+par.mpc_eps)-c[t,:])/par.mpc_eps

        
        # if t< par.T-1:
        #     if t+1 > par.Tr: #after pension
        #         m[t+1,:] = par.R*a[t,:]/(par.G*par.L[t])+1
        #         p[t+1,:] = np.log(par.G)+np.log(par.L[t])+p[t,:]
        #         y[t+1,:] = p[t+1,:]
        #     else:            #before pension
        #         m[t+1,:] = par.R*a[t,:]/(par.G*par.L[t]*sim.psi[t+1,:])+sim.xi[t+1,:]
        #         p[t+1,:] = np.log(par.G)+np.log(par.L[t])+p[t,:]+np.log(sim.psi[t+1,:])
        #         y[t+1,:] = p[t+1,:]+np.log(sim.xi[t+1,:])
              
        