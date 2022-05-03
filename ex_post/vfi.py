import numpy as np
from numba import njit, prange
from scipy import optimize

 # consav
from consav import linear_interp # for linear interpolation
from consav import golden_section_search # for optimization in 1D
import tools
import utility

# a. define objective function
def obj_bellman(sol, par, c,t,m):
        """ value of choice of c used in vfi """

        # a. end-of-period assets
        a = m-c

        # b. next-period cash-on-hand
        still_working_next_period = t+1 <= par.TR-1
        if still_working_next_period:
            fac = par.G*par.L[t]*par.psi_vec
            w = par.w
            xi = par.xi_vec
        else:
            fac = par.G*par.L[t]
            w = 1
            xi = 1

        m_plus = (par.R/fac)*a + xi            

        # c. continuation value
        if still_working_next_period:
            inv_v_plus = np.zeros(m_plus.size)
            tools.interp_linear_1d(sol.m[t+1,:],sol.inv_v[t+1,:],m_plus,inv_v_plus)
        else:
            inv_v_plus = tools.interp_linear_1d_scalar(sol.m[t+1,:],sol.inv_v[t+1,:],m_plus)
        v_plus = 1/inv_v_plus
        
        # d. value-of-choice
        total = utility(c,par) + par.beta*np.sum(w*fac**(1-par.rho)*v_plus)
        return -total

# b. solve bellman equation        
def VFI(t,sol,par):
    """solve bellman equation using vfi"""

    # a. last period (= consume all)
    sol.m[-1,:] = par.grid_m[-1,:]
    sol.c[-1,:] = sol.m[-1,:]
    for i,c in enumerate(sol.c[-1,:]):
        sol.inv_v[-1,i] = 1.0/utility.func(c,par)
    
    # b. before last period
    for t in reversed(range(par.T-1)):
        for i_m in range(par.Nm):

            m = par.grid_m[t,i_m]

            obj = lambda c: obj_bellman(sol,par,c,t,m)
            result = optimize.minimize_scalar(obj,method='bounded',bounds=(0,m))

            sol.c[t,i_m] = result.x
            sol.inv_v[t,i_m]= -1.0/result.fun
        
        # save grid for m
        sol.m[t,:] = par.grid_m[t,:]