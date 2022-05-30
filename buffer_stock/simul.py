import numpy as np
from numba import njit,prange

from consav import linear_interp
import tools

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
    
    for t in range(par.T):
        
        grid_m = sol_m[t,:]      
        grid_c = sol_c[t,:]
        
        c[t,:] = tools.interp_linear_1d(grid_m,grid_c,m[t,:])
        a[t,:] = m[t,:] - c[t,:]
        
        if t< par.T-1:
            if t+1 > par.Tr: #after pension
                m[t+1,:] = par.R*a[t,:]/(par.G*par.L[t])+1
                p[t+1,:] = np.log(par.G)+np.log(par.L[t])+p[t,:]
                y[t+1,:] = p[t+1,:]
            else:       #before pension
                m[t+1,:] = par.R*a[t,:]/(par.G*par.L[t]*sim.psi[t+1,:])+sim.xi[t+1,:]
                p[t+1,:] = np.log(par.G)+np.log(par.L[t])+p[t,:]+np.log(sim.psi[t+1,:])
                y[t+1,:] = p[t+1,:]+np.log(sim.xi[t+1,:])
              
        mpc[t,:] = (tools.interp_linear_1d(grid_m,grid_c,m[t,:]+par.mpc_eps)-c[t,:])/par.mpc_eps
        #if t == 4:
            #print((tools.interp_linear_1d(grid_m,grid_c,4+par.mpc_eps)-tools.interp_linear_1d(grid_m,grid_c,4))/par.mpc_eps)
    