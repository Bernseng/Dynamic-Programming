import numpy as np
from numba import njit, prange

from consav import linear_interp

import trans

@njit(parallel=True)
def life_cycle(par,sol,sim):
    """ simulate model with parallization over households """

    # unpack (helps numba)
    m = sim.m
    p = sim.p
    y = sim.y
    c = sim.c
    a = sim.a

    sol_c = sol.c
    sol_m = sol.m

    # loop over first households and then time
    for i in prange(par.simN):
        for t in range(par.T):
            
            # a. solution
            grid_m = sol_m[t,:]
            grid_c = sol_c[t,:]
            
            # b. consumption
            c[i,t] = linear_interp.interp_1d(grid_m,grid_c,m[i,t])
            a[i,t] = m[i,t] - c[i,t]

            # c. next-period
            if t == 0:
                m[i,t] = sim.m0[i]
                p[i,t] = sim.p0[i]
            else:
                m[i,t] = trans.m_plus_func(a[i,t-1],sim.xi[i,t],sim.psi[i,t],par,t-1)
                p[i,t] = trans.p_plus_func(p[i,t-1],sim.psi[i,t],par,t-1)
            
            y[i,t] = trans.y_plus_func(p[i,t],sim.xi[i,t],par,t-1)
            
            
            