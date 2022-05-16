# -*- coding: utf-8 -*-
"""
Prints welfare figure.

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

import matplotlib.pyplot as plt
                     
plt.rcParams.update({'font.size': 8,
                     'font.family': 'STIXGeneral',
                     'mathtext.fontset': 'stix'})

Color1 = (3.0/255.0,103.0/255.0,166.0/255.0)
Color2 = (242.0/255.0,62.0/255.0,46.0/255.0)
Color3 = (3.0/255.0,166.0/255.0,166.0/255.0)
Color4 = (242.0/255.0,131.0/255.0,68.0/255.0)
Color5 = (242.0/255.0,100.0/255.0,48.0/255.0)
ColorGrey = (65.0/255.0,68.0/255.0,81.0/255.0)                     
 
from plot_sim import SimStruct

def fig_explanation():
    
    path = os.getcwd() + '\\Welfare'
    if not os.path.exists(path):
        os.makedirs(path)
        
    # a. load baseline data    
    par = pickle.load(open('baseline_med_0\\Data\\par.p', 'rb' ))                       
    filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/sim')    
    base_sim = SimStruct()
         
    for x in ['pop_age', 'pop_c', 'pop_P']:
                  
        exec('temp = g[\'' + x + '\']')         
        exec('base_sim.' + x + ' = np.copy(temp[...])')  
                   
    f.close()

    # b. load high lambda data    
    par = pickle.load(open('lamb_5\\Data\\par.p', 'rb' ))                       
    filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/sim')    
    alt_sim = SimStruct()
         
    for x in ['pop_age', 'pop_c', 'pop_P']:
                  
        exec('temp = g[\'' + x + '\']')         
        exec('alt_sim.' + x + ' = np.copy(temp[...])')  
                   
    f.close()

    # c. calculate stuff
    mean_c = np.zeros((120,2))
    std_c = np.zeros((120,2))   
    for t in xrange(120):
        
        I = base_sim.pop_age == t
        mean_c[t,0] = np.mean(base_sim.pop_c[I])
        std_c[t,0] = np.std(base_sim.pop_c[I])
        
        I = alt_sim.pop_age == t
        mean_c[t,1] = np.mean(alt_sim.pop_c[I])
        std_c[t,1] = np.std(alt_sim.pop_c[I])
        
    
    # d. mean figure
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   

    ax = fig.add_subplot(1,1,1)        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    ax.plot(np.arange(120), mean_c[:,0],
            ls='none', lw=2, 
            marker='o', markersize=6, 
            color=Color1,     
            markerfacecolor = Color1, markeredgecolor = Color1,             
            label='$\\lambda = 0.03$')

    ax.plot(np.arange(120), mean_c[:,1],
            ls='none', lw=2, 
            marker='o', markersize=6, 
            color=Color2,     
            markerfacecolor = Color2, markeredgecolor = Color2,             
            label='$\\lambda = 0.99$')
             
    legend = ax.legend(loc='lower right', fontsize=10, numpoints=1)
    
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)
             
    fig.tight_layout()    
    fig.savefig(os.getcwd() + '//Welfare/welfare_mean_c.pdf')
    plt.close('all')      
   
    # e. std figure
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
   
    ax = fig.add_subplot(1,1,1)        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    ax.plot(np.arange(120), std_c[:,0],
            ls='none', lw=2, 
            marker='o', markersize=6, 
            color=Color1,     
            markerfacecolor = Color1, markeredgecolor = Color1,             
            label='$\\lambda = 0.03$')

    ax.plot(np.arange(120), std_c[:,1],
            ls='none', lw=2, 
            marker='o', markersize=6, 
            color=Color2,     
            markerfacecolor = Color2, markeredgecolor = Color2,             
            label='$\\lambda = 0.99$')
             
    legend = ax.legend(loc='lower right', fontsize=10, numpoints=1)
    
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)
             
    fig.tight_layout()    
    fig.savefig(os.getcwd() + '//Welfare/welfare_std_c.pdf')
    plt.close('all')      
            
    return
    
def fig_welfare(parlamb, vals, vals_alt, names_alt):
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    
    ax = fig.add_subplot(1,1,1)        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    ax.plot(vals[:, 0], vals[:, 1], 
             ls='-', lw=2, 
             marker='o', markersize=4, 
             color=Color1,     
             markerfacecolor = Color1, markeredgecolor = Color1,             
             label=r'$\tau$')
    ax.set_xlabel('$\lambda$', fontsize=10)
    ax.set_ylabel('percent', color=Color1, fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.set_xlim([0, 1])
    
    plt.axvline(x=parlamb, ymin=0, ymax=100, linewidth=3, color='black', alpha=0.5)
    
    ax.plot(parlamb, vals_alt[0], 
             ls='none', lw=2, 
             marker='o', markersize=6, 
             color=Color4,     
             markerfacecolor = Color4, markeredgecolor = Color4,             
             label=names_alt[0])

    ax.plot(parlamb, vals_alt[1], 
             ls='none', lw=2, 
             marker='^', markersize=6, 
             color=Color4,     
             markerfacecolor = Color4, markeredgecolor = Color4,         
             label=names_alt[1])        
         
    legend = ax.legend(loc='lower right', fontsize=10, numpoints=1)
    
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)
    
    fig.tight_layout()    
    fig.savefig(os.getcwd() + '//Welfare/welfare.pdf')
    plt.close('all')      
    
def find_vals(par):

    # a. load data                      
    filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/sim')    
    sim = SimStruct()
         
    for x in ['pop_age', 'pop_c', 'pop_P']:
                  
        exec('temp = g[\'' + x + '\']')         
        exec('sim.' + x + ' = np.copy(temp[...])')  
                   
    f.close()
    
    I = sim.pop_c[:] == 0
    print(' = 0: ' + str(I.sum()))
    I = sim.pop_c[:] < 0
    print(' < 0: ' + str(I.sum()))
    
    I = sim.pop_c[:] > 0
    beta = (par.beta+par.pi_death)**sim.pop_age[I]    
    C = sim.pop_c[I]*sim.pop_P[I]*(par.Gamma**sim.pop_age[I])    
    u = C**(1.0-par.rho)/(1.0-par.rho)
    u_all = beta*u
    u_exp = np.mean(u_all)        
    
    I = u[:] >= 0
    print(' >= 0: ' + str(I.sum()))
    
    return np.array([u_exp])

def allfigs():
    
    path = os.getcwd() + '\\Welfare'
    if not os.path.exists(path):
        os.makedirs(path)
        
    files = glob.glob(path + '//*.pdf')
    for f in files:
        os.remove(f)
        
    fig_explanation()
        
    versions = [['lamb','']]
      
    print('base')          
    par = pickle.load(open('baseline_med_0\\Data\\par.p', 'rb' ))            
    base_vals = find_vals(par)

    names_alt = []
    vals_alt = np.zeros(2)

    print('alt 0')
    par_alt = pickle.load(open('sigma_xi_3\\Data\\par.p', 'rb' ))            
    alt_vals = find_vals(par_alt)
    vals_alt[0] = ((base_vals[0] / alt_vals[0])**(1.0/(1.0-par.rho))-1.0)*100
    names_alt.append(r'$\sigma_{\xi}$: 0.20 $\rightarrow$ 0.30 ')
        
    print('alt 1')
    par_alt = pickle.load(open('u_ast_3\\Data\\par.p', 'rb' ))             
    alt_vals = find_vals(par_alt)
    vals_alt[1] = ((base_vals[0] / alt_vals[0])**(1.0/(1.0-par.rho))-1.0)*100
    names_alt.append(r'$u_{\ast}$: 0.07 $\rightarrow$ 0.14') 
        
    print(vals_alt)
    for version in versions:

        name = version[0] 
            
        print(name)          
    
        folders = glob.glob(os.getcwd() + '//' + name + '*')
        folders.append('baseline_med_0')
            
        vals = np.zeros((len(folders),2))
        
        for i, folder in enumerate(folders):
            
            print((i))
            par = pickle.load(open( folder + '\\Data\\par.p', 'rb' ))      
            
            exec('vals[i, 0] = par.' + name)
                     
            vals[i, 1:] = find_vals(par)
            vals[i, 1] = ((base_vals[0] / vals[i, 1])**(1.0/(1.0-par.rho))-1.0)*100
   
        parlamb = vals[-1, 0] 
                       
        I = np.argsort(vals[:, 0])
        vals[:, 0] = vals[I, 0]
        vals[:, 1] = vals[I, 1]
        
        fig_welfare(parlamb, vals, vals_alt, names_alt)
             
    print(vals_alt)
    print(vals)
            
    return

# main #
if __name__ == '__main__':
        
    allfigs()