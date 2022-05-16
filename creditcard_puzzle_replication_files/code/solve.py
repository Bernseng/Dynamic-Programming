# -*- coding: utf-8 -*-
"""
Solution module.

Version: 1.0.
@author: Jeppe Druedahl, 2017.
"""

from __future__ import division, print_function, absolute_import
import os
clear = lambda: os.system('cls')
clear()

import time
import glob
import h5py
import numpy as np
                    
from cffi import FFI
ffi = FFI()

# compile c #
os.system("gcc -fopenmp -o3 -std=gnu99 -c cfuncs\\vfi.c -o cfuncs\\vfi.o")
os.system("gcc -fopenmp -shared -o cfuncs\\vfi.so cfuncs\\vfi.o")

ffi.cdef(r''' void vfi(double *v_T, double *d_T, double *Delta_T, double *c_T, double *n_T, double *a_T, double *vt_T, double *ucon_c_d_T,
            double *db, double *nb_T, double *kappa, int *i_nb_f_T, int *Nnb_db_T,
            double *psi, double *xit, double *w, double *w_scale,
            double beta, double rho, double Gamma,
            double ra, double rd, double lambda, double eta, double varphi, double chi,
            double n_max, double phi_n, double eps,
            int T, int tcon, int Nu, int Nx, int Ndb, int Nnb_max, int Npsi, int  Nxit, int Nn, long int Nd, double d_stepsize); ''')
    
lib = ffi.dlopen("cfuncs\\vfi.so")
    
# functions #    
def vfi(p):    

    print("solve.py")
    
    # a. allocate memory
    t1 = time.time()
    
    d = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, p.Nnb_max))
    Delta = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, p.Nnb_max))

    c = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, p.Nnb_max))
    n = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, p.Nnb_max))
    a = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, p.Nnb_max))

    v = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, p.Nnb_max))
    vt = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, p.Nnb_max))
    
    ucon_c_d = np.zeros((2, p.Nu, p.Nx, p.Nd, p.Nnb_max)) # due to memory

    # b. create pointers
    p_d = ffi.cast('double *', d.ctypes.data)            
    p_Delta = ffi.cast('double *', Delta.ctypes.data)                
    p_c = ffi.cast('double *', c.ctypes.data)            
    p_a = ffi.cast('double *', a.ctypes.data)            
    p_n = ffi.cast('double *', n.ctypes.data)            
    p_v = ffi.cast('double *', v.ctypes.data)                
    p_vt = ffi.cast('double *', vt.ctypes.data)            

    p_ucon_c_d = ffi.cast('double *', ucon_c_d.ctypes.data)              

    p_db = ffi.cast('double *', p.db_vec.ctypes.data)        
    p_nb = ffi.cast('double *', p.nb_mat.ctypes.data)     
    p_i_nb_f = ffi.cast('int *', p.i_nb_f.ctypes.data)             
    p_Nnb_db = ffi.cast('int *', p.Nnb_db.ctypes.data)       

    p_kappa = ffi.cast('double *', p.kappa.ctypes.data)        
    
    p_psi = ffi.cast('double *', p.psi_vec.ctypes.data)
    p_xit = ffi.cast('double *', p.xit_vec.ctypes.data)       
    p_w = ffi.cast('double *', p.w.ctypes.data)  
    p_w_scale = ffi.cast('double *', p.w_scale.ctypes.data) 

    print(" - memory allocated (time = {:3.1f})".format(time.time()-t1))
    
    # c. run    
    t1 = time.time()
        
    lib.vfi(p_v, p_d, p_Delta, p_c, p_n, p_a, p_vt, p_ucon_c_d,
            p_db, p_nb, p_kappa, p_i_nb_f, p_Nnb_db,
            p_psi, p_xit, p_w, p_w_scale,
            p.beta, p.rho, p.Gamma,
            p.ra, p.rd, p.lamb, p.eta, p.varphi, p.chi,
            p.n_max, p.phi_n, p.epsilon,
            p.T, p.tcon, p.Nu, p.Nx, p.Ndb, p.Nnb_max, p.Npsi, p.Nxit, p.Nn, np.int32(p.Nd), p.d_stepsize)
              
    print(" - model solved (time = {:3.1f})".format(time.time()-t1))
    
    # d. save results
    t1 = time.time()
    
    files = glob.glob(os.getcwd() + '\\' + p.name + '\\Data\\*_' + str(p.het) + '.h5')
    for f in files:
        os.remove(f)
        
    filename = os.getcwd() + '\\' + p.name + '\\Data\\Data_' + str(p.het) + '.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/sol')
    
    for x in ['v', 'vt', 'd', 'Delta', 'c', 'a', 'n', 'ucon_c_d']:

        if x in g:
            del g[x]
            
        g[x] = eval('' + x + '[0, :, :, :, :]')
        
    f.close()

    print(" - data saved (time = {:3.1f})".format(time.time()-t1))

    return
              
# main #
if __name__ == "__main__":
    pass    