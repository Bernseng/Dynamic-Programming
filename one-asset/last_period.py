import numpy as np
from numba import njit, prange

from consav import golden_section_search
from consav.grids import nonlinspace_jit

import utility

@njit
def last_period(sol,par):
    
    #unpack
    m = sol.m
    c = sol.c
    inv_v = sol.inv_v
    
    # last period consume all       
    m[-1,:] = nonlinspace_jit(0,par.a_max,par.Na+1,par.m_phi)
    c[-1,:] = sol.m[-1,:]
    inv_v[-1,0] = 0
    inv_v[-1,1:] = 1.0/utility.func(sol.c[-1,1:],par)