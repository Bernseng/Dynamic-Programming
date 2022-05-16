# -*- coding: utf-8 -*-
"""
Prints solution figures.

Version: 1.0.
@author: Jeppe Druedahl, 2017.
"""

# imports #
from __future__ import division, print_function, absolute_import
import os
clear = lambda: os.system('cls')
clear()

import glob
import h5py
import cPickle as pickle
import numpy as np
import time

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

from itertools import product

def policy_and_value_figs(par, i_u, i_x):

    # a. load data        
    filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
    t = 0
        
    f = h5py.File(filename, 'a')
    g = f.require_group('/sol')
    vt = np.copy(g['vt'])
    vt = vt[i_u, i_x, :, :]
    
    # b. which figures to print?
    names = ['vt', 'd', 'Delta', 'c', 'a', 'n']
    names_latex = [r'$\tilde{v}_t$', r'$d_t$', r'$\Delta_t$', r'$c_t$', r'$a_t$', r'$n_t$']

    assert len(names) == len(names_latex)

    # c. which d-indexes?
    i_db_vec = np.array([0, 30, 40, 50, 55, par.Ndb-1])
    num_i_db = i_db_vec.size

    # d. figures loop
    nb = par.nb_mat[t, i_u, i_x, :, :]

    for j, name in enumerate(names):
        
        exec("y = np.copy(g[\'" + name + "\'])")
        
        fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
        ax = fig.add_subplot(1,1,1)
        
        ax.spines["top"].set_linewidth(0.5) 
        ax.spines["bottom"].set_linewidth(0.5) 
        ax.spines["right"].set_linewidth(0.5)  
        ax.spines["left"].set_linewidth(0.5)  
        ax.tick_params(which='both', color=ColorGrey)
        ax.grid(True, zorder=0, color=ColorGrey)
    
        for i in xrange(num_i_db):

            i_db = i_db_vec[i]

            if names[j] in ['d', 'a']:
                i_f = par.i_nb_f[t, i_u, i_x, i_db]+1
            else:
                i_f = par.i_nb_f[t, i_u, i_x, i_db]

            # not including all elements in nb_end
            i_l = par.i_nb_l[t, i_u, i_x, i_db] - (len(par.nb_end)-1)

            nb_now = nb[i_db, i_f:i_l]
            y_now = y[i_u, i_x, i_db, i_f:i_l]
            vt_now = vt[i_db, i_f:i_l]
            
            if names[j] in ['d', 'Delta', 'c', 'a', 'n']:
                I = vt_now > 0
                nb_now = nb_now[I]
                y_now = y_now[I]
                
            ax.plot(nb_now, y_now,
                    lw=1, ls='-', marker='o', markersize=2.5, markeredgecolor='none',
                    color=cm.jet(1.*i/num_i_db),
                    label=r"$\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")

        ax.set_xlabel(r'$\bar{n}_t$', fontsize=10)
        ax.set_ylabel(names_latex[j], fontsize=10)

        if names[j] in ['c']:                
            legend = ax.legend(loc='lower right', fontsize=8)
        elif names[j] in ['d']:                
            legend = ax.legend(loc='upper right', fontsize=8)
        else:                
            legend = ax.legend(loc='best', fontsize=8)            
        
        frame = legend.get_frame()        
        frame.set_linewidth(0.5)   
        frame.set_edgecolor(ColorGrey)
        
        [line.set_zorder(3) for line in ax.lines]
        fig.tight_layout()
        fig.savefig(os.getcwd() + '\\' + par.name + "\\Graphs\\Sol\\" + names[j] + "_u" + str(i_u) + "_x" + str(i_x) + ".pdf")
        
        plt.close('all')
        
    # c. special figure with c and d
    c = np.copy(g['c'])
    d = np.copy(g['d'])
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    for i in [0, num_i_db-1]:

        i_db = i_db_vec[i]

        i_f = par.i_nb_f[t, i_u, i_x, i_db]
        i_l = par.i_nb_l[t, i_u, i_x, i_db] - (len(par.nb_end)-1)

        nb_now = nb[i_db, i_f:i_l]
        c_now = c[i_u, i_x, i_db, i_f:i_l]       
        vt_now = vt[i_db, i_f:i_l]            
        I = vt_now > 0
        nb_c = nb_now[I]
        c_now = c_now[I]

        nb_now = nb[i_db, i_f+1:i_l]
        d_now = d[i_u, i_x, i_db, i_f+1:i_l]
        vt_now = vt[i_db, i_f+1:i_l]            
        I = vt_now > 0
        nb_d = nb_now[I]
        d_now = d_now[I]
        
        if i == 0: 

            ax.plot(nb_c, c_now,
                    lw=1, ls='-', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color1,
                    label=r"$c_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")
    
            ax.plot(nb_d, d_now,
                    lw=1, ls='-', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color2,
                    label=r"$d_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")
        else:
            
            ax.plot(nb_c, c_now,
                    lw=1, ls='--', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color1,
                    label=r"$c_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")
    
            ax.plot(nb_d, d_now,
                    lw=1, ls='--', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color2,
                    label=r"$d_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")            
                        
        ax.set_xlabel(r'$\bar{n}_t$', fontsize=10)
        ax.set_ylabel(r'$d_t$, $c_t$', fontsize=10)

        legend = ax.legend(loc='best', fontsize=8)            
        
        frame = legend.get_frame()        
        frame.set_linewidth(0.5)   
        frame.set_edgecolor(ColorGrey)
        
        [line.set_zorder(3) for line in ax.lines]
        fig.tight_layout()
        fig.savefig(os.getcwd() + '\\' + par.name + "\\Graphs\\Sol\\cd_u" + str(i_u) + "_x" + str(i_x) + ".pdf")
        
        plt.close('all')

    # c. special figure with a and d
    a = np.copy(g['a'])
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    for i in [0, num_i_db-1]:

        i_db = i_db_vec[i]

        i_f = par.i_nb_f[t, i_u, i_x, i_db]
        i_l = par.i_nb_l[t, i_u, i_x, i_db] - (len(par.nb_end)-1)

        nb_now = nb[i_db, i_f:i_l]
        a_now = a[i_u, i_x, i_db, i_f:i_l]       
        vt_now = vt[i_db, i_f:i_l]            
        I = vt_now > 0
        nb_a = nb_now[I]
        a_now = a_now[I]

        nb_now = nb[i_db, i_f+1:i_l]
        d_now = d[i_u, i_x, i_db, i_f+1:i_l]
        vt_now = vt[i_db, i_f+1:i_l]            
        I = vt_now > 0
        nb_d = nb_now[I]
        d_now = d_now[I]
        
        if i == 0: 

            ax.plot(nb_a, a_now,
                    lw=1, ls='-', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color1,
                    label=r"$a_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")
    
            ax.plot(nb_d, d_now,
                    lw=1, ls='-', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color2,
                    label=r"$d_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")
        else:
            
            ax.plot(nb_a, a_now,
                    lw=1, ls='--', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color1,
                    label=r"$c_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")
    
            ax.plot(nb_d, d_now,
                    lw=1, ls='--', marker='', markersize=2.5, markeredgecolor='none',
                    color=Color2,
                    label=r"$d_t$, $\bar{d}_t = " + str(round(par.db_vec[i_db]*10)/10) + "$")            
                        
        ax.set_xlabel(r'$\bar{n}_t$', fontsize=10)
        ax.set_ylabel(r'$d_t$, $a_t$', fontsize=10)

        legend = ax.legend(loc='best', fontsize=8)            
        
        frame = legend.get_frame()        
        frame.set_linewidth(0.5)   
        frame.set_edgecolor(ColorGrey)
        
        [line.set_zorder(3) for line in ax.lines]
        fig.tight_layout()
        fig.savefig(os.getcwd() + '\\' + par.name + "\\Graphs\\Sol\\ca_u" + str(i_u) + "_x" + str(i_x) + ".pdf")
        
        plt.close('all')
        
    f.close()     
    return
    
    
def puzzle_figs(p, i_u, i_x):
        
    filename = os.getcwd() + '\\' + p.name + '\\Data\\Data_' + str(p.het) + '.h5'
        
    f = h5py.File(filename, 'a')
    g = f.require_group('/sol')

    vt = np.copy(g['vt'])    
    d = np.copy(g['d'])
    a = np.copy(g['a'])
    n = np.copy(g['n'])
    
    f.close()

    vt = np.ravel(vt[i_u, i_x, :, :])
    d = np.ravel(d[i_u, i_x, :, :])
    a = np.ravel(a[i_u, i_x, :, :])
    n = np.ravel(n[i_u, i_x, :, :])
    
    db_mat = np.zeros((p.Ndb, p.Nnb_max))
    for i_db in xrange(p.Ndb):
        db_mat[i_db, :] = p.db_vec[i_db]
    
    db = np.ravel(db_mat)    
    nb = np.ravel(p.nb_mat[0, i_u, i_x, :, :])
    
    Ivt = vt > 0.0
    Id = d > 0.0
    Ia = a > 0.0

    Ipuzzle = (Ivt == True) & (Id == True) & (Ia == True) 
    Isaver = (Ivt == True) & (Ipuzzle == False) & (Ia == True) & (Id == False) 
    Iborrower = (Ivt == True) & (Ipuzzle == False) & (Ia == False) & (Id == True)
    Icorner = (Ivt == True) & (Ipuzzle == False) & (Ia == False) & (Id == False)
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
    
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    puzzle = ax.scatter(db[Ipuzzle], nb[Ipuzzle], c=Color1, edgecolor='none', s=10, marker='.')       
    borrower = ax.scatter(db[Iborrower], nb[Iborrower], c=Color2, edgecolor='none', s=8, marker='v')    
    saver = ax.scatter(db[Isaver], nb[Isaver], c=Color3, edgecolor='none', s=8, marker='^')
    corner = ax.scatter(db[Icorner], nb[Icorner], c=Color4, edgecolor='none', s=6, marker='s')
    
    legend = ax.legend((puzzle, borrower, saver, corner),
                       ('$a_t > 0$, $d_t > 0$', '$a_t = 0$, $d_t > 0$', '$a_t > 0$, $d_t = 0$','$a_t = 0$, $d_t = 0$'),
                       scatterpoints=1,
                       loc='lower left',
                       ncol=2,
                       fontsize=8,
                       markerscale=3)

    frame = legend.get_frame()        
    frame.set_linewidth(0.5)   
    frame.set_edgecolor(ColorGrey)
         
    ax.set_xlabel(r'$\bar{d}_t$', fontsize=10)
    ax.set_ylabel(r'$\bar{n}_t$', fontsize=10)
                    
    ax.set_xlim([-0.1, 2 + 0.05])            
    ax.set_ylim([p.kappa[0, i_u, i_x, -1]-0.5, 2 + 0.1])
    
    fig.tight_layout()    
    fig.savefig(os.getcwd() + '\\' + p.name + "\\Graphs\\Sol\\groups_u" + str(i_u) + "_x" + str(i_x) + ".pdf")    
    plt.close('all')

    return
        
########    
# MENU #
########
          
def allfigs(par):

    files = glob.glob(os.getcwd() + '\\' + par.name + '//Graphs/Sol/*.pdf')
    for f in files:
        os.remove(f) 
            
    print("plot_sol.py")

    t1 = time.time()
    for (i_u, i_x) in product(xrange(par.Nu), xrange(par.Nx)):
        puzzle_figs(par, i_u, i_x)
    
    print(" - groups ({:3.1f})".format(time.time()-t1))  
        
    t1 = time.time()
    for (i_u, i_x) in product(xrange(par.Nu), xrange(par.Nx)):
        policy_and_value_figs(par, i_u, i_x)
    
    print(" - value and policy functions ({:3.1f})".format(time.time()-t1))
            
    return
    
# main #
if __name__ == '__main__':
    pass