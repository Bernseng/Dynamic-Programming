import numpy as np

import tools
import trans

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
