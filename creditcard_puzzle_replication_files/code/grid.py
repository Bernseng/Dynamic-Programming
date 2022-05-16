# -*- coding: utf-8 -*-
"""
Creates the grid.

0. Generic Functions
1. Shock Functions
2. State Space Functions
3. Tabs and Figs
4. Full Grid

Version: 1.0.
@author: Jeppe Druedahl, 2017.
"""

from __future__ import division, print_function, absolute_import
import os
clear = lambda: os.system('cls')
clear()

import time
import glob
import cPickle as pickle
import numpy as np
from itertools import product

from scipy.optimize import minimize_scalar

import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams.update({'font.size': 8,
                     'font.family': 'STIXGeneral',
                     'mathtext.fontset': 'stix'})

Color1 = (3.0/255.0,103.0/255.0,166.0/255.0)
Color2 = (242.0/255.0,62.0/255.0,46.0/255.0)
Color3 = (3.0/255.0,166.0/255.0,166.0/255.0)
Color4 = (242.0/255.0,131.0/255.0,68.0/255.0)
Color5 = (242.0/255.0,100.0/255.0,48.0/255.0)
ColorGrey = (65.0/255.0,68.0/255.0,81.0/255.0)

import settings

from cffi import FFI
ffi = FFI()

# compile c #
os.system("gcc -o3 -std=gnu99 -c cfuncs\\check_bounds_func.c -o cfuncs\\check_bounds_func.o")
os.system("gcc -shared -o cfuncs\\check_bounds_func.so cfuncs\\check_bounds_func.o")

ffi.cdef(r''' void check_bounds_func(int *check, int t, int i_u, int i_x, double db, double nb, double *db_plus,
                                     double *kappa_plus,
                                     double *psi, double *xit,
                                     double Gamma, double ra, double rd, double lambda,
                                     double eta, double varphi,
                                     int Nu, int Nx, int Ndb, long int Nd, int Npsi, int Nxit, double d_stepsize); ''')

lib = ffi.dlopen("cfuncs\\check_bounds_func.so")


########################
# 0. GENERIC FUNCTIONS #
########################

def nonlinspace(x_min, x_max, n, phi, x_end=[]):
    """
    Returns a vector of unequally distributed points between x_min and x_max
    (both included) with the points in x_end inserted at the end.

    INPUT::

     x_min:     minimum value
     x_max:     maximum value (x_max > x_min)
     n:         number of points to return >= 1
     phi:       curvature parameter >= 1
     x_end:     (optional) values inserted after the end of x

    OUTPUT::

     y:         asc sorted 1D nd.array, len=n

    A higher phi *increases* the frequency of points close to x_min.
    For phi = 1 the stepsize is constant.

    Examples:

    >>> y = nonlinspace(1, 9, 3, 1)
    >>> y
    array([ 1.,  5.,  9.])
    >>> y = nonlinspace(1, 9, 3, 2)
    >>> y
    array([ 1.,  3.,  9.])
    """

    assert x_max > x_min
    assert n >= 2
    assert phi >= 1

    # 1. n's
    n_extra = len(x_end)
    n_base = n-n_extra

    # 2. recursion
    y = np.empty(n)

    y[0] = x_min
    for i in xrange(1, n_base):
        y[i] = y[i-1] + (x_max-y[i-1]) / (n_base-i)**phi

    # 3. add at end
    for i in xrange(0, n_extra):
        y[n_base + i] = x_end[i]

    assert np.all(np.diff(y) > 0)

    return y


######################
# 1. SHOCK FUNCTIONS #
######################

def gauss_hermite(sigma, mu, num):
    """
    Returns nodes and weights for a log-guassian shock
    using Gauss-Hermite quadrature.

    INPUT::

     sigma:     std., >= 0
     mu:        mean
     num:       number of nodes, == 1 or even > 0

    OUTPUT::

     nodes      1D nd.array, num
     weights    1D nd.array, num

    The weights sum to one.

    Example:

    >>> nodes, weights = gauss_hermite(1, 0, 8)
    >>> np.allclose(sum(weights), 1)
    True
    """

    assert sigma >= 0
    assert num >= 1
    assert num == 1 or num % 2 == 0, "num cannot be uneven"

    if num == 1 or sigma == 0:
        nodes = np.array([1])
        weights = np.array([1])
    else:
        nodes, weights = np.polynomial.hermite.hermgauss(num)
        nodes *= 2**0.5
        nodes = np.exp(sigma*nodes + mu)
        weights *= np.pi**-0.5

    return nodes, weights


def create_shocks(p):
    """
    Returns node and weight *vectors* for two *uncorelated*
    log-gaussian shocks using double Gauss-Hermite quadrature
    adding a low probability *extreme* value to the second shock.

    OUTPUT::

        psi:          desc sorted, psi nodes, 1D nd.array, Npsi*(Nxi+1)
        xit:          asc sorted, xit nodes, 1D nd.array, Npsi*(Nxi+1)
        w:            associated weights, 1D nd.array, len=Npsi*(Nxi+1)
        w:            scaled weights, 1D nd.array, len=Npsi*(Nxi+1)
        
    Calls: gauss_hermite()
    """

    ####################
    # 1. income shocks #
    ####################

    # gauss hermite
    psi, psi_w = gauss_hermite(p.sigma_psi, p.mu_psi, p.Npsi)
    xi, xi_w = gauss_hermite(p.sigma_xi, p.mu_xi, p.Nxi)

    psi = np.copy(psi[::-1])
    psi_w = np.copy(psi_w[::-1])

    # xi_tilde
    Nxit = p.Nxi+p.Nu-1
    xit = np.empty(Nxit)
    xit_w = np.empty((p.Nu, Nxit)) 
        
    if p.Nu == 2:
        
        xit[0] = p.mu
        xit[1:] = (xi-p.mu*p.u_ast) / (1-p.u_ast)
    
        xit_w[0, 0] = p.pi_uw
        xit_w[0, 1:] = xi_w*(1-p.pi_uw)
        xit_w[1, 0] = p.pi_uu
        xit_w[1, 1:] = xi_w*(1-p.pi_uu)
    
    else:
        
        xit[:] = xi
        xit_w[0, :] = xi_w
    
    # check sorting
    
    if p.Npsi > 1:
        assert np.all(np.diff(psi) < 0), psi

    if p.Nxi > 1:
        assert np.all(np.diff(xit[1:]) > 0), xit


    ####################
    # 2. credit shocks #
    ####################

    p.pi_x_lose = np.empty((p.Nu))
    p.pi_x_gain = np.empty((p.Nu))
        
    if p.x_additive:
        
        p.pi_x_gain = np.empty((p.Nu))
        p.pi_x_gain[0] = p.pi_x_gain_tot
        p.pi_x_gain[1] = p.pi_x_gain[0] + p.chi_gain
        
        p.pi_x_lose = np.empty((p.Nu))    
        p.pi_x_lose[0] = p.pi_x_lose_tot
        p.pi_x_lose[1] = p.pi_x_lose[0] + p.chi_lose
        
    else:
        
        def f(x, u_ast, fac, target):
            return (target - ((1-u_ast)*x + u_ast * fac * x))**2

        res = minimize_scalar(f, args=(p.u_ast, p.chi_lose, p.pi_x_lose_tot),
                              bounds=(0, 1), method='bounded')
        
        p.pi_x_lose[0] = res.x          
        p.pi_x_lose[1] = p.chi_lose * res.x
        
        # gain
        p.pi_x_gain[0] = p.pi_x_gain_tot
        p.pi_x_gain[1] = p.pi_x_gain_tot
     
    # transition matrix
    pi_x_mat = np.empty((p.Nu,p.Nx,p.Nx))       
    for i_u in xrange(p.Nu):        
        pi_x_mat[i_u, 0, 0]  = 1.0 - p.pi_x_lose[i_u]
        pi_x_mat[i_u, 0, 1]  = p.pi_x_lose[i_u]
        pi_x_mat[i_u, 1, 1]  = 1.0 - p.pi_x_gain[i_u]
        pi_x_mat[i_u, 1, 0]  = p.pi_x_gain[i_u]
        
        
    # 4. total weight vector
    w = np.empty((p.Nu,p.Nx,p.Nx,p.Npsi,Nxit))
    w_scale = np.empty((p.Nu,p.Nx,p.Nx,p.Npsi,Nxit))
    
    for (i_u, i_x, i_x_plus, i_psi, i_xit) in product(xrange(p.Nu), 
                                                      xrange(p.Nx), 
                                                      xrange(p.Nx), 
                                                      xrange(p.Npsi), 
                                                      xrange(Nxit)):

        if i_xit == 0 and p.Nu > 1:
            i_u_plus = 1
        else:
            i_u_plus = 0
            
        pi_x = pi_x_mat[i_u_plus, i_x, i_x_plus]   
        
        w_now =  pi_x*psi_w[i_psi]*xit_w[i_u, i_xit]
        scale_now = (p.Gamma*psi[i_psi])**(1.0-p.rho)
        
        w[i_u, i_x, i_x_plus, i_psi, i_xit] = w_now
        w_scale[i_u, i_x, i_x_plus, i_psi, i_xit] =  w_now*scale_now

    for i_u, i_x in product(xrange(p.Nu), xrange(p.Nx)):
        testsum = np.sum(w[i_u, i_x, :, :, :])
        assert np.allclose(testsum, 1.0), testsum

    return psi, xit, psi_w, xit_w, w, w_scale, Nxit


############################
# 2. STATE SPACE FUNCTIONS #
############################

def create_db(p):
    
    db1 = nonlinspace(0, 2, p.Ndb-20, 1.1)
    db2 = np.linspace(2+0.1, p.Upsilon, 20)
    db = np.unique(np.hstack((db1, db2)))    
    
    assert db[0] == 0
    assert db[-1] == p.Upsilon
    assert len(db) == p.Ndb

    return db


def create_kappa(p):
    
    t1 = time.time()

    kappa = np.empty((p.T, p.Nu, p.Nx, p.Ndb))
    check = np.empty(1, dtype=int)

    # 1. last period
    for (i_u, i_x) in product(xrange(p.Nu), xrange(p.Nx)):
        kappa[p.T-1, i_u, i_x, :] = 0
        kappa[p.T-1, i_u, i_x, :] = 0
    
    # 2. remaining periods (backwards)
    p_db_plus = ffi.cast('double *', p.db_vec.ctypes.data)
    p_psi     = ffi.cast('double *', p.psi_vec.ctypes.data)
    p_xit     = ffi.cast('double *', p.xit_vec.ctypes.data)
    p_check   = ffi.cast('int *', check.ctypes.data)

    t = p.T-2
    while(t >= 0):
                
        for (i_u, i_x, i_db) in product(xrange(p.Nu), xrange(p.Nx), xrange(p.Ndb-1, -1, -1)):
            
            t1 = time.time()
    
            kappa_plus = kappa[t+1, :, :, :]
                
            p_kappa_plus = ffi.cast('double *', kappa_plus.ctypes.data)
                    
            db = p.db_vec[i_db]
    
            # a. liquidity is directly bounding
            if i_x == 1:
                varphi_now = 0
            elif i_u == 1:
                varphi_now = p.chi*p.varphi
            else:
                varphi_now = p.varphi  
    
            #nb_low = -np.fmax(db, varphi_now/2)
            nb_low = -np.fmax(db, varphi_now/(1+p.eta))            
            nb_high = 0
    
            # b. restrict it to be a decreasing function
            if i_db < p.Ndb-1:
                nb_low = np.fmax(nb_low, kappa[t, i_u, i_x, i_db+1])
    
            # c. inner recursion
            lib.check_bounds_func(p_check, t, i_u, i_x, db, nb_low, p_db_plus,
                                  p_kappa_plus, 
                                  p_psi, p_xit,
                                  p.Gamma, p.ra, p.rd, p.lamb, 
                                  p.eta, varphi_now,
                                  p.Nu, p.Nx, p.Ndb, np.int32(p.Nd), p.Npsi, p.Nxit, p.d_stepsize)
    
            if check[0] == 1:
                
                kappa[t, i_u, i_x, i_db] = nb_low
                
            else:
                    
                if t == 0:
                    pass #print("INNER")
                    
                while(True):
    
                    nb_mid = (nb_high+nb_low)/2
    
                    lib.check_bounds_func(p_check, t, i_u, i_x, db, nb_mid, p_db_plus,
                                          p_kappa_plus, 
                                          p_psi, p_xit,
                                          p.Gamma, p.ra, p.rd, p.lamb, 
                                          p.eta, varphi_now,
                                          p.Nu, p.Nx, p.Ndb, np.int32(p.Nd), p.Npsi, p.Nxit, p.d_stepsize)

                    if check[0] == 1:
                        nb_high = nb_mid
                    else:
                        nb_low = nb_mid
    
                    if nb_high-nb_low < p.epsilon:
                        kappa[t, i_u, i_x, i_db] = nb_high
                        break
                    
        # decrement
        t -= 1
            
    return kappa


def create_nb(p):

    Nnb_max = p.Ndb + p.Nnb
    nb = np.zeros((p.T, p.Nu, p.Nx, p.Ndb, Nnb_max))
    i_nb_f = np.empty((p.T, p.Nu, p.Nx, p.Ndb), dtype=int)
    i_nb_l = np.empty((p.T, p.Nu, p.Nx, p.Ndb), dtype=int)
    Nnb_db = np.empty((p.T, p.Nu, p.Nx, p.Ndb), dtype=int)

    for (i_u, i_x, t) in product(xrange(p.Nu), xrange(p.Nx), xrange(p.T)):

        # a. unique minimum values
        uniq_kappa = np.unique(p.kappa[t, i_u, i_x, :])
    
        # b. common nodes
        nb_nodes = nonlinspace(np.max(uniq_kappa), p.nb_max,
                               p.Nnb, p.phi_nb, p.nb_end)
    
        # c. stacking and sorting
        nb_full = np.unique(np.hstack((uniq_kappa, nb_nodes)))
        nb_full.sort()
        assert np.all(np.diff(nb_full) >= 0)
    
        # create grid
        Nnb_full = len(nb_full)
        for i in xrange(p.Ndb):
    
            Nnb_db[t, i_u, i_x, i] = (nb_full >= p.kappa[t,i_u, i_x, i]).sum()
            i_nb_f[t, i_u, i_x, i] = Nnb_full - Nnb_db[t, i_u, i_x, i]
            i_nb_l[t, i_u, i_x, i] = Nnb_full
        
            nb[t, i_u, i_x, i, i_nb_f[t, i_u, i_x, i]:i_nb_l[t, i_u, i_x, i]] = nb_full[i_nb_f[t, i_u, i_x, i]:i_nb_l[t, i_u, i_x, i]]
        
            assert np.all(np.diff(nb[t, i_u, i_x, i, i_nb_f[t, i_u, i_x, i]:i_nb_l[t, i_u, i_x, i]]) >= 0)
            assert nb[t, i_u, i_x, i, i_nb_f[t, i_u, i_x, i]] == p.kappa[t, i_u, i_x, i]

        assert np.all(np.diff(Nnb_db[t, i_u, i_x, :]) >= 0)

    return nb, Nnb_db, i_nb_f, i_nb_l, Nnb_max


####################
# 3. TABS AND FIGS #
####################

def print_vec(var, name):
    print("\n\nVariable: {}\n".format(name))
    for i in range(len(var)):
        print(" {}[{}] = {}".format(name, i, var[i]))
    return


def fig_state_space_t(t, i_u, i_x, p):

    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)

    for i in range(p.Ndb):
        nb_vec = p.nb_mat[t, i_u, i_x, i, p.i_nb_f[t, i_u, i_x, i]:p.i_nb_l[t, i_u, i_x, i]]
        ax.plot(p.db_vec[i]*np.ones(p.Nnb_db[t, i_u, i_x, i]),
                nb_vec, 'o', color=Color1, markeredgecolor='none', markersize=3)

    ax.set_xlabel(r'$\bar{d}_t$', fontsize=10)
    ax.set_ylabel(r'$\bar{n}_t$', fontsize=10)

    ax.set_xlim([0, p.db_vec[-1]+0.5])
    nb_vec = np.ravel(p.nb_mat[t, i_u, i_x, :, :])
    ymin = np.min(nb_vec)  
    ax.set_ylim([ymin-0.5, p.nb_end[-1]+0.5])

    [line.set_zorder(3) for line in ax.lines]
    fig.tight_layout()
    version = '_t' + str(np.int(t)) + '_u' + str(i_u) + '_x' + str(i_x)
    fig.savefig(os.getcwd() + '\\' + p.name + '\\Graphs\\Grid\\state_space' + version + '.pdf')
    
    plt.close('all')
    return
    
def fig_state_space_kappa(i_u, i_x, p):

    if p.T > 24:
        t_vec = p.T - np.array([1, 4, 8, 12, 16, p.T-8, p.T-4, p.T])
    else:
        t_vec = p.T - np.floor(nonlinspace(1, p.T, 5, 1))

    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)

    for counter, t in enumerate(t_vec):

        ax.plot(p.db_vec, p.kappa[t, i_u, i_x, :],
                lw=1, ls='-', marker='.', markersize=3,
                color=cm.jet(1.*counter/len(t_vec)),
                label=r'$t = ' + str(p.T-np.int(t)) + '$')

        ax.set_xlabel(r'$\bar{d}_t$', fontsize=10)
        ax.set_ylabel(r'$\bar{n}_t$', fontsize=10)
        ax.set_xlim([0, p.db_vec[-1]])
        ax.set_ylim([p.kappa[0, i_u, i_x, -1]-0.5, 0.1])

    legend = ax.legend(loc='best', fontsize=8)
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)

    [line.set_zorder(3) for line in ax.lines]
    fig.tight_layout()
    version = '_u' + str(i_u) + '_x' + str(i_x)
    fig.savefig(os.getcwd() + '\\' + p.name  + '\\Graphs\\Grid\\state_space_kappa' + version + '.pdf')

    plt.close('all')
    return

def allfigs(p):
    
    t1 = time.time()
    
    for (i_u, i_x) in product(xrange(p.Nu), xrange(p.Nx)):

        fig_state_space_kappa(i_u, i_x, p)
        fig_state_space_t(0, i_u, i_x, p)

    print(" - figures (time = {:3.1f})".format(time.time()-t1))
    return


################
# 4. FULL GRID #
################

def create(p):
    
    if not os.path.exists(os.getcwd() + '\\' + p.name):
        os.makedirs(os.getcwd() + '\\' + p.name )
        os.makedirs(os.getcwd() + '\\' + p.name  + '\\Graphs')
        os.makedirs(os.getcwd() + '\\' + p.name  + '\\Data')
        os.makedirs(os.getcwd() + '\\' + p.name  + '\\Graphs\\Grid') 
        os.makedirs(os.getcwd() + '\\' + p.name  + '\\Graphs\\Sol')
        os.makedirs(os.getcwd() + '\\' + p.name  + '\\Graphs\\Sim')
                
    files = glob.glob(os.getcwd() + '\\' + p.name + '//Graphs//Grid//*')
    for f in files:
        os.remove(f)

    with open("log_check_bounds.txt", "w") as text_file:
        pass
    
    print("grid.py")
    
    # a. build grid
    t1 = time.time()

    p.mu_psi = -0.5 * p.sigma_psi**2
    p.mu_xi = -0.5 * p.sigma_xi**2       
    p.pi_uw = (1-p.pi_uu)/(1-p.u_ast)*p.u_ast
    p.rd = p.ra + p.drd
    p.n_max = p.nb_max 
        
    p.Nd = int(np.floor(p.Upsilon/p.d_stepsize))
    p.psi_vec, p.xit_vec, p.psi_w, p.xit_w, p.w, p.w_scale, p.Nxit = create_shocks(p)
    p.db_vec = create_db(p)
    p.kappa = create_kappa(p)
    p.nb_mat, p.Nnb_db, p.i_nb_f, p.i_nb_l, p.Nnb_max = create_nb(p)

    tcon = p.T-2
    while(tcon >= 0 and np.max(np.abs(p.kappa[tcon+1, :, :, :] - p.kappa[tcon, :, :, :])) != 0):
        tcon -= 1
    p.tcon = tcon

    print(" - built (time = {:3.1f})".format(time.time()-t1))
    
    # b. grid figures
    #allfigs(p)
    
    # c. save par
    t1 = time.time()
    pickle.dump(p, open( os.getcwd() + '\\' + p.name  + '\\Data\\par.p', 'wb' ))
    print(" - data saved (time = {:3.1f})".format(time.time()-t1))
            
    return
    
# main    
if __name__ == "__main__":
    pass