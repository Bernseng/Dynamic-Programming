import numpy as np
from numba import njit, prange

 # consav
from consav import linear_interp # for linear interpolation

# local modules
import trans
import utility

@njit(parallel=True)
def lifecycle(sim,sol,par):
    """ simulate full life-cycle """

    # unpack
    p = sim.p
    y = sim.y
    n = sim.n
    m = sim.m
    c = sim.c
    d = sim.d
    a = sim.a
    c_bump = sim.c_bump
    mpc = sim.mpc
    discrete = sim.discrete
    
    for t in range(par.T):
        for i in prange(par.simN):
            
            # a. beginning of period states
            if t == 0:
                p[t,i] = sim.p0[i]
                n[t,i] = sim.d0[i]
                m[t,i] = sim.a0[i]
            else:
                p[t,i] = trans.p_plus_func(p[t-1,i],sim.psi[t,i],par,t-1)
                n[t,i] = trans.n_plus_func(d[t-1,i],par,sim.z[t])
                m[t,i] = trans.m_plus_func(a[t-1,i],p[t,i],sim.xi[t,i],par,t)
            
            y[t,i] = p[t,i] * sim.xi[t,i]

            # b. optimal choices and post decision states

            optimal_choice(i,t,p[t,i],n[t,i],m[t,i],discrete[t,i:],d[t,i:],c[t,i:],a[t,i:],sol,par,mpc[t,i:],c_bump[t,i:])
            
@njit            
def optimal_choice(i,t,p,n,m,discrete,d,c,a,sol,par,mpc,c_bump):

    x = trans.x_plus_func(m,n,par)
    x_mpc = trans.x_plus_func(m+par.mpc_eps,n,par)

    # a. discrete choice
    inv_v_keep = linear_interp.interp_3d(par.grid_p,par.grid_n,par.grid_m,sol.inv_v_keep[t],p,n,m)    
    inv_v_keep_mpc = linear_interp.interp_3d(par.grid_p,par.grid_n,par.grid_m,sol.inv_v_keep[t],p,n,m+par.mpc_eps)
    inv_v_adj = linear_interp.interp_2d(par.grid_p,par.grid_x,sol.inv_v_adj[t],p,x)    
    inv_v_adj_mpc = linear_interp.interp_2d(par.grid_p,par.grid_x,sol.inv_v_adj[t],p,x_mpc)    
    adjust = inv_v_adj > inv_v_keep
    adjust_mpc = inv_v_adj_mpc > inv_v_keep_mpc

    
    # b. continuous choices
    if adjust:

        discrete[0] = 1
        
        d[0] = linear_interp.interp_2d(
            par.grid_p,par.grid_x,sol.d_adj[t],
            p,x)

        c[0] = linear_interp.interp_2d(
            par.grid_p,par.grid_x,sol.c_adj[t],
            p,x)
        
        tot = d[0]+c[0]

        if tot > x: 
            d[0] *= x/tot
            c[0] *= x/tot
            a[0] = 0.0
        else:
            a[0] = x - tot

        # calculate mpc
        if par.cross_compute:

            if adjust==adjust_mpc:
                c_bump[0] = linear_interp.interp_2d(par.grid_p,par.grid_x,sol.c_adj[t],
                    p,x_mpc)
            else:
                c_bump[0] = linear_interp.interp_3d(
                par.grid_p,par.grid_n,par.grid_m,sol.c_keep[t],
                p,n,m+par.mpc_eps)
        else:
            c_bump[0] = linear_interp.interp_2d(par.grid_p,par.grid_x,sol.c_adj[t],
                p,x_mpc)

        tot = d[0]+c_bump[0]

        if tot > x+par.mpc_eps: 
            d[0] *= (x+par.mpc_eps)/tot
            c_bump[0] *= (x+par.mpc_eps)/tot
            a[0] = 0.0
        else:
            a[0] = (x+par.mpc_eps) - tot

        mpc[0] = (c_bump[0] - c[0]) / par.mpc_eps

    else: 
            
        discrete[0] = 0

        d[0] = n

        c[0] = linear_interp.interp_3d(
            par.grid_p,par.grid_n,par.grid_m,sol.c_keep[t],
            p,n,m)

        if c[0] > m: 
            c[0] = m
            a[0] = 0.0
        else:
            a[0] = m - c[0]
        
        # calculate mpc
        if par.cross_compute:
            if adjust==adjust_mpc:
                c_bump[0] = linear_interp.interp_3d(
                    par.grid_p,par.grid_n,par.grid_m,sol.c_keep[t],
                    p,n,m+par.mpc_eps)
            else:
                c_bump[0] = linear_interp.interp_2d(par.grid_p,par.grid_x,sol.c_adj[t],
                    p,x_mpc)
        else: 
            c_bump[0] = linear_interp.interp_3d(
                par.grid_p,par.grid_n,par.grid_m,sol.c_keep[t],
                p,n,m+par.mpc_eps)

        if c_bump[0] > (m+par.mpc_eps): 
            c_bump[0] = m+par.mpc_eps
            a[0] = 0.0
        else:
            a[0] = m + par.mpc_eps - c_bump[0]

        mpc[0] = (c_bump[0] - c[0]) / par.mpc_eps



@njit            
def euler_errors(sim,sol,par):

    # unpack
    euler_error = sim.euler_error
    euler_error_c = sim.euler_error_c
    
    for i in prange(par.simN):
        
        discrete_plus = np.zeros(1)
        d_plus = np.zeros(1)
        c_plus = np.zeros(1)        
        c_bump_plus = np.zeros(1)
        a_plus = np.zeros(1)

        for t in range(par.T-1):

            constrained = sim.a[t,i] < par.euler_cutoff
            
            if constrained:

                euler_error[t,i] = np.nan
                euler_error_c[t,i] = np.nan
                continue

            else:

                RHS = 0.0
                for ishock in range(par.Nshocks):
                        
                    # i. shocks
                    psi = par.psi[ishock]
                    psi_w = par.psi_w[ishock]
                    xi = par.xi[ishock]
                    xi_w = par.xi_w[ishock]

                    # ii. next-period states
                    p_plus = trans.p_plus_func(sim.p[t,i],psi,par)
                    n_plus = trans.n_plus_func(sim.d[t,i],par)
                    m_plus = trans.m_plus_func(sim.a[t,i],p_plus,xi,par,t)

                    # iii. weight
                    weight = psi_w*xi_w

                    # iv. next-period choices

                    optimal_choice(t+1,p_plus,n_plus,m_plus,discrete_plus,d_plus,c_plus,a_plus,sol,par,c_bump_plus)

                    # v. next-period marginal utility

                    RHS += weight*par.beta*par.R*utility.marg_func(c_plus[0],d_plus[0],par)
                

                euler_error[t,i] = sim.c[t,i] - utility.inv_marg_func(RHS,sim.d[t,i],par)

                euler_error_c[t,i] = sim.c[t,i]

@njit(parallel=True)
def calc_utility(sim,sol,par):
    """ calculate utility for each individual """

    # unpack
    u = sim.utility
    
    for t in range(par.T):
        for i in prange(par.simN):
            
            u[i] += par.beta**t*utility.func(sim.c[t,i],sim.d[t,i],par)
            
