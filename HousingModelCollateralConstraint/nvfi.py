import numpy as np
from numba import njit, prange

# consav
from consav import linear_interp # for linear interpolation
from consav import golden_section_search

# local modules
import utility

########
# keep #
########

@njit
def obj_keep(c,n,b,inv_w,grid_a,d_ubar,alpha,rho,theta,omega):
    """ evaluate bellman equation """

    # a. end-of-period assets
    a = b-c-omega*n
    
    # b. continuation value
    w = -1.0/linear_interp.interp_1d(grid_a,inv_w,a)

    # c. total value
    value_of_choice = utility.func_nopar(c,n,d_ubar,alpha,rho,theta) + w

    return -value_of_choice # we are minimizing


@njit(parallel=True)
def solve_keep(t,sol,par):
    """solve bellman equation for keepers using nvfi"""

    # unpack output
    inv_v = sol.inv_v_keep[t]
    inv_marg_u = sol.inv_marg_u_keep[t]
    c = sol.c_keep[t]

    # unpack input
    inv_w = sol.inv_w[t]
    grid_a = par.grid_a
    d_ubar = par.d_ubar
    alpha = par.alpha
    rho = par.rho
    theta = par.theta
    omega = par.omega[t]

    # loop over outer states
    for i_p in prange(par.Np):
        for i_n in range(par.Nn):
            
            # outer states
            n = par.grid_n[i_n]

            # loop over b state
            for i_b in range(par.Nb):
                
                # a. cash-on-hand
                b = par.grid_b[i_b]
                if i_b == 0:
                    c[i_p,i_n,i_b] = 0
                    inv_v[i_p,i_n,i_b] = 0
                    if par.do_marg_u:
                        inv_marg_u[i_p,i_n,i_b] = 0        
                    continue

                # b. optimal choice
                c_low = np.fmin(b/2,1e-8)
                c_high = b
                c[i_p,i_n,i_b] = golden_section_search.optimizer(obj_keep,c_low,c_high,
                                                                          args=(n,b,inv_w[i_p,i_n],grid_a[i_n,:],
                                                                                d_ubar,alpha,rho,theta,omega),
                                                                          tol=par.tol)

                # c. optimal value
                v = -obj_keep(c[i_p,i_n,i_b],n,b,inv_w[i_p,i_n],grid_a[i_n,:],d_ubar,alpha,rho,theta,omega)
                inv_v[i_p,i_n,i_b] = -1/v
                if par.do_marg_u:
                    inv_marg_u[i_p,i_n,i_b] = 1/utility.marg_func_nopar(c[i_p,i_n,i_b],n,d_ubar,alpha,rho,theta)

#######
# adj #
#######

@njit
def obj_adj(d,x,inv_v_keep,grid_n,grid_b,omega):
    """ evaluate bellman equation """

    # a. cash-on-hand
    b = x-d-omega*d

    # b. durables
    n = d
    
    # c. value-of-choice
    return -linear_interp.interp_2d(grid_n,grid_b,inv_v_keep,n,b)  # we are minimizing

@njit(parallel=True)
def solve_adj(t,sol,par):
    """solve bellman equation for adjusters using nvfi"""

    # unpack output
    inv_v = sol.inv_v_adj[t]
    inv_marg_u = sol.inv_marg_u_adj[t]
    d = sol.d_adj[t]
    c = sol.c_adj[t]

    # unpack input
    inv_v_keep = sol.inv_v_keep[t]
    c_keep = sol.c_keep[t]
    grid_n = par.grid_n
    grid_m = par.grid_m
    grid_b = par.grid_b
    d_ubar = par.d_ubar
    alpha = par.alpha
    rho = par.rho
    theta = par.theta
    omega = par.omega[t]

    # loop over outer states
    for i_p in prange(par.Np):
            
        # loop over x state
        for i_x in range(par.Nx):
            
            # a. cash-on-hand
            x = par.grid_x[i_x]
            if i_x == 0:
                d[i_p,i_x] = 0
                c[i_p,i_x] = 0
                inv_v[i_p,i_x] = 0
                if par.do_marg_u:
                    inv_marg_u[i_p,i_x] = 0        
                continue

            # b. optimal choice
            d_low = np.fmin(x/2,1e-8)
            d_high = np.fmin(x/(1-omega),par.n_max)
            d[i_p,i_x] = golden_section_search.optimizer(obj_adj,d_low,d_high,args=(x,inv_v_keep[i_p],grid_n,grid_b,omega),tol=par.tol)

            # c. optimal value
            b = x - d[i_p,i_x] + omega*d[i_p,i_x]
            c[i_p,i_x] = linear_interp.interp_2d(par.grid_n,par.grid_b,c_keep[i_p],d[i_p,i_x],b)
            inv_v[i_p,i_x] = -obj_adj(d[i_p,i_x],x,inv_v_keep[i_p],grid_n,grid_b,omega)
            if par.do_marg_u:
                inv_marg_u[i_p,i_x] = 1/utility.marg_func_nopar(c[i_p,i_x],d[i_p,i_x],d_ubar,alpha,rho,theta)