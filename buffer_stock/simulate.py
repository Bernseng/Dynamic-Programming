import numpy as np
from numba import njit, prange

from consav import linear_interp

import trans

#@njit(parallel=True)
def life_cycle(par,sol,sim):
    """ simulate model with parallization over households """

    # unpack (helps numba)
    m = sim.m
    p = sim.p
    y = sim.y
    c = sim.c
    a = sim.a
    mpc = sim.mpc

    sol_c = sol.c
    sol_m = sol.m

    # loop over first households and then time
    for t in prange(par.T):
        for i in range(par.simN):
            
            # a. solution
            grid_m = sol_m[t,:]
            grid_c = sol_c[t,:]
            
            # b. consumption
            c[t,i] = linear_interp.interp_1d(grid_m,grid_c,m[t,i])
            a[t,i] = m[t,i] - c[t,i]

            # c. next-period
            if t == 0:
                m[t,i] = sim.m0[i]
                p[t,i] = sim.p0[i]
            else:
                m[t,i] = trans.m_plus_func(a[t-1,i],sim.xi[t,i],sim.psi[t,i],par,t-1)
                p[t,i] = trans.p_plus_func(p[t-1,i],sim.psi[t,i],par,t-1)
            
            y[t,i] = trans.y_plus_func(p[t,i],sim.xi[t,i],par,t-1)
            
            mpc[t,i] = (linear_interp.interp_1d(grid_m,grid_c,m[t,i]+par.mpc_eps)-c[t,i])/par.mpc_eps
            
            
            