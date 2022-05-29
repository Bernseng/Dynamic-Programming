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
    
    sol_c = sol.c
    sol_m = sol.m 
    
    sim.m[0,:] = sim.m0
    sim.p[0,:] = sim.p0
    
    for t in range(par.T):
        
        grid_m = sol_m[t,:]      
        grid_c = sol_c[t,:]
        
        c[t,:] = tools.interp_linear_1d(grid_m,grid_c,m[t,:])
        a[t,:] = m[t,:] - c[t,:]
        
        if t< par.T-1:
            if t+1 > par.Tr: #after pension
                sim.m[t+1,:] = par.R*sim.a[t,:]/(par.G*par.L[t])+1
                sim.p[t+1,:] = np.log(par.G)+np.log(par.L[t])+sim.p[t,:]
                sim.y[t+1,:] = sim.p[t+1,:]
            else:       #before pension
                sim.m[t+1,:] = par.R*sim.a[t,:]/(par.G*par.L[t]*sim.psi[t+1,:])+sim.xi[t+1,:]
                sim.p[t+1,:] = np.log(par.G)+np.log(par.L[t])+sim.p[t,:]+np.log(sim.psi[t+1,:])
                sim.y[t+1,:] = sim.p[t+1,:]+np.log(sim.xi[t+1,:])
              
            
    
    
    