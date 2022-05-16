# -*- coding: utf-8 -*-
"""
Simulation module.

version: 1.0.
@author: Jeppe Druedahl, 2017.
"""

# imports #
from __future__ import division, print_function, absolute_import
import os
clear = lambda: os.system('cls')
clear()

import time
import h5py
import cPickle as pickle
import numpy as np
import plot_sim
                     
from cffi import FFI
ffi = FFI()

# compile c #
os.system("gcc -o3 -std=gnu99 -c cfuncs\\sim_func.c -o cfuncs\\sim_func.o")
os.system("gcc -shared -o cfuncs\\sim_func.so cfuncs\\sim_func.o")

ffi.cdef(r''' void sim_func(int N, int T,
         int *pop_age, int *pop_u, int *pop_x, double *pop_db, double *pop_nb,
         double *pop_d, double *pop_c,
         double *pop_n, double *pop_a,
         double *pop_death, double *pop_u_val, double *pop_x_val,
         double *pop_psi, double *pop_xi,
         double *pop_P, double *pop_Y, double *pop_util,
         int Nu, int Nx, int Ndb, double *db,
         int Nnb_max, double *nb_T, int *i_nb_f_T, int *Nnb_db_T,
         long int Nd, double d_stepsize,
         double *d_T, double *ucon_c_d_T,
         double rho, double Gamma, double mu,
         double ra, double rd, double lambda, double eta, double varphi, double chi,
         double pi_uw, double pi_uu, double *pi_x_lose, double *pi_x_gain, double pi_death); ''')
                                   
lib = ffi.dlopen("cfuncs\\sim_func.so")
        
# functions #
def simulate(p):

    print("sim.py")
    
    # a. load policy functions    
    t1 = time.time()
    
    filename = os.getcwd() + '\\' + p.name + '\\Data\\Data_' + str(p.het) + '.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/sol')  

    d = np.copy(g['d'])
    ucon_c_d = np.copy(g['ucon_c_d'])    
    
    f.close()

    p_d_T = ffi.cast('double *', d.ctypes.data)
    p_ucon_c_d_T = ffi.cast('double *', ucon_c_d.ctypes.data)
    
    print(" - policy functions loaded (time = {:3.1f})".format(time.time()-t1))
        
    # b. allocate memory
    t1 = time.time()

    pop_age = np.zeros((p.simT+1, p.N), dtype=int)      
    pop_u = np.zeros((p.simT+1, p.N), dtype=int)  
    pop_x = np.zeros((p.simT+1, p.N), dtype=int)  
    
    pop_db = np.zeros((p.simT+1, p.N))    
    pop_nb = np.zeros((p.simT+1, p.N))
    pop_d = np.zeros((p.simT, p.N))
    pop_c = np.zeros((p.simT, p.N))
    pop_n = np.zeros((p.simT, p.N))
    pop_a = np.zeros((p.simT, p.N))

    pop_P = np.zeros((p.simT+1, p.N))
    pop_Y = np.zeros((p.simT+1, p.N))
    
    pop_util = np.zeros((p.simT, p.N))
    
    print(" - memory allocated (time = {:3.1f})".format(time.time()-t1))
    
    # c. draw shocks
    t1 = time.time()
    
    np.random.seed(seed=1986)

    # psi
    if p.Npsi > 1:
        
        temp = np.exp(np.random.normal(p.mu_psi, p.sigma_psi, p.N*(p.simT+1)))

        I = temp < np.min(p.psi_vec)
        temp[I] = np.min(p.psi_vec)
        I = temp > np.max(p.psi_vec)
        temp[I] = np.max(p.psi_vec)
        
        pop_psi = temp.reshape((p.simT+1, p.N))

    else:
        
        pop_psi = np.ones((p.simT+1, p.N))    
        
    # xi
    if p.Nxi > 1:

        temp = np.exp(np.random.normal(p.mu_xi, p.sigma_xi, p.N*(p.simT+1)))
        temp = (temp-p.mu*p.u_ast) / (1.0-p.u_ast)
        
        I = temp < np.min(p.xit_vec[1:])
        temp[I] = np.min(p.xit_vec[1:])
        I = temp > np.max(p.xit_vec[1:])
        temp[I] = np.max(p.xit_vec[1:])
        
        pop_xi = temp.reshape((p.simT+1, p.N))    

    else:
        
        pop_xi = np.ones((p.simT+1, p.N))  
        
    print(" - random normal numbers (time = {:3.1f})".format(time.time()-t1))

    t1 = time.time()
        
    filename = 'random.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/random')
    
    if 'pop_death' not in g:
             
        # death
        pop_death = np.zeros((p.simT, p.N))
        for t in xrange(p.simT):
            pop_death[t, :] = np.linspace(0.0, 1.0, p.N)
            np.random.shuffle(pop_death[t, :])
            
        # u_val
        pop_u_val = np.zeros((p.simT+1, p.N))
        for t in xrange(p.simT+1):
            pop_u_val[t, :] = np.linspace(0.0, 1.0, p.N)
            np.random.shuffle(pop_u_val[t, :])
            
        # x_val    
        pop_x_val = np.zeros((p.simT+1, p.N))        
        for t in xrange(p.simT+1):
            pop_x_val[t, :] = np.linspace(0.0, 1.0, p.N)
            np.random.shuffle(pop_x_val[t, :])
            
        g['pop_death'] = pop_death
        g['pop_u_val'] = pop_u_val
        g['pop_x_val'] = pop_x_val
        
    else:

        pop_death = np.copy(g['pop_death'])
        pop_u_val = np.copy(g['pop_u_val'])
        pop_x_val = np.copy(g['pop_x_val'])
        
    f.close()

    print(" - random uniform numbers (time = {:3.1f})".format(time.time()-t1))
       
    # d. pointers
    if p.no_x == 1:
        p.pi_xu = 0.0
            
    t1 = time.time()
            
    p_pop_age = ffi.cast('int *', pop_age.ctypes.data)        
    p_pop_u = ffi.cast('int *', pop_u.ctypes.data)        
    p_pop_x = ffi.cast('int *', pop_x.ctypes.data)        
        
    p_pop_db = ffi.cast('double *', pop_db.ctypes.data)        
    p_pop_nb = ffi.cast('double *', pop_nb.ctypes.data)        
    p_pop_d = ffi.cast('double *', pop_d.ctypes.data)        
    p_pop_c = ffi.cast('double *', pop_c.ctypes.data)        
    p_pop_n = ffi.cast('double *', pop_n.ctypes.data)        
    p_pop_a = ffi.cast('double *', pop_a.ctypes.data)        

    p_pop_P = ffi.cast('double *', pop_P.ctypes.data)        
    p_pop_Y = ffi.cast('double *', pop_Y.ctypes.data)

    p_pop_util = ffi.cast('double *', pop_util.ctypes.data)
    
    p_pop_death = ffi.cast('double *', pop_death.ctypes.data)
    p_pop_u_val = ffi.cast('double *', pop_u_val.ctypes.data)  
    p_pop_x_val = ffi.cast('double *', pop_x_val.ctypes.data)      
    p_pop_psi = ffi.cast('double *', pop_psi.ctypes.data)        
    p_pop_xi = ffi.cast('double *', pop_xi.ctypes.data)        

    p_db = ffi.cast('double *', p.db_vec.ctypes.data)
    p_nb_T = ffi.cast('double *', p.nb_mat.ctypes.data)
    p_i_nb_f_T = ffi.cast('int *', p.i_nb_f.ctypes.data)
    p_Nnb_db_T = ffi.cast('int *', p.Nnb_db.ctypes.data)

    p_pi_x_lose = ffi.cast('double *', p.pi_x_lose.ctypes.data)
    p_pi_x_gain = ffi.cast('double *', p.pi_x_gain.ctypes.data)
    
    print(" - pointers set (time = {:3.1f})".format(time.time()-t1))
     
    # e. simulate
    t1 = time.time()
    
    lib.sim_func(p.N, p.simT,
                 p_pop_age, p_pop_u, p_pop_x, p_pop_db, p_pop_nb,
                 p_pop_d, p_pop_c,
                 p_pop_n, p_pop_a,
                 p_pop_death, p_pop_u_val, p_pop_x_val,
                 p_pop_psi, p_pop_xi,
                 p_pop_P, p_pop_Y, p_pop_util,
                 p.Nu, p.Nx, p.Ndb, p_db,
                 p.Nnb_max, p_nb_T, p_i_nb_f_T, p_Nnb_db_T,
                 np.int32(p.Nd), p.d_stepsize,
                 p_d_T, p_ucon_c_d_T,
                 p.rho, p.Gamma, p.mu,
                 p.ra, p.rd, p.lamb, p.eta, p.varphi, p.chi,
                 p.pi_uw, p.pi_uu, p_pi_x_lose, p_pi_x_gain, p.pi_death)              
    
    print(" - simulation done (time = {:3.1f})".format(time.time()-t1))

    # f. save   
    t1 = time.time()
     
    filename = os.getcwd() + '\\' + p.name  + '\\Data\\Data_' + str(p.het) + '.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/sim')    
    
    for x in ['pop_age', 'pop_u', 'pop_x', 'pop_db', 'pop_nb', 
              'pop_d', 'pop_c', 'pop_n', 'pop_a', 
              'pop_death', 'pop_psi', 'pop_xi',
              'pop_P', 'pop_Y', 'pop_util']:

        if x in g:
            del g[x]
            
        if x in ['pop_db', 'pop_nb', 
                 'pop_d', 'pop_c', 'pop_n', 'pop_a',
                 'pop_death', 'pop_psi', 'pop_xi',
                 'pop_P', 'pop_Y', 'pop_util']:
        
            g[x] = eval('np.single(' + x + '[p.simBurnIn:p.simT])')
        
        else:
            
            g[x] = eval('' + x + '[p.simBurnIn:p.simT]')

    f.close()  

    print(" - data saved (time = {:3.1f})".format(time.time()-t1))


    return
     
     
# main #
if __name__ == '__main__':
    pass