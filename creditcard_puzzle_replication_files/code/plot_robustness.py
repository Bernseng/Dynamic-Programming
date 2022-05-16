# -*- coding: utf-8 -*-
"""
Prints robustness figures.

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
                     
plt.rcParams.update({'font.size': 9,
                     'font.family': 'STIXGeneral',
                     'mathtext.fontset': 'stix'})

Color1 = (3.0/255.0,103.0/255.0,166.0/255.0)
Color2 = (242.0/255.0,62.0/255.0,46.0/255.0)
Color3 = (3.0/255.0,166.0/255.0,166.0/255.0)
Color4 = (242.0/255.0,131.0/255.0,68.0/255.0)
Color5 = (242.0/255.0,100.0/255.0,48.0/255.0)
ColorGrey = (65.0/255.0,68.0/255.0,81.0/255.0)                     
 
from plot_sim import SimStruct
from plot_sim import calc_sim
    
def fig_robustness(parbase, name, name_latex, vals, legendpars):
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   

    ax1 = fig.add_subplot(1,1,1)        
    ax1.spines["top"].set_linewidth(0.5) 
    ax1.spines["bottom"].set_linewidth(0.5) 
    ax1.spines["right"].set_linewidth(0.5)  
    ax1.spines["left"].set_linewidth(0.5)  
    ax1.tick_params(which='both', color=ColorGrey)
    ax1.grid(True, zorder=0, color=ColorGrey)
    
    ax1.plot(vals[:, 0], vals[:, 1], 
             ls='-', lw=2, 
             marker='o', markersize=4, 
             color=Color1,     
             markerfacecolor = Color1, markeredgecolor = Color1,             
             label=r'puzzle share')

    ax1.plot(vals[:, 0], vals[:, 2], 
             ls='--', lw=2, 
             marker='o', markersize=4, 
             color=Color1,     
             markerfacecolor = Color1, markeredgecolor = Color1,             
             label=r'borrower share')

             
    ax1.set_xlabel(name_latex, fontsize=11)
    ax1.set_ylabel('percent', color=Color1, fontsize=11)
    ax1.tick_params(axis='both', which='major', labelsize=11)    
    ax1.set_ylim([0, 60])
    xmin = vals[0, 0] - 0.1 * (vals[-1, 0] - vals[0, 0])
    xmax = vals[-1, 0] + 0.1 * (vals[-1, 0] - vals[0, 0])
    ax1.set_xlim([xmin, xmax])
    
    plt.axvline(x=parbase, ymin=0, ymax=100, linewidth=3, color='black', alpha=0.5)

    for tl in ax1.get_yticklabels():
        tl.set_color(Color1)

    ax2 = ax1.twinx()
           
    ax2.plot(vals[:, 0], vals[:, 3], 
             ls='-.', lw=2, 
             marker='o', markersize=4,
             color=Color2,   
             markerfacecolor = Color2, markeredgecolor = Color2,
             label=r'credit card debt, $D_t$') 

    ax2.plot(vals[:, 0], vals[:, 4], 
             ls=':', lw=2, 
             marker='o', markersize=4, 
             color=Color2,   
             markerfacecolor = Color2, markeredgecolor = Color2,
             label=r'liquid assets, $A_t$') 
             
    ax2.set_ylabel(r'$D_t$, $A_t$', color=Color2, fontsize=11)
    ax2.tick_params(axis='both', which='major', labelsize=11)
    ax2.set_ylim([-0.1, 0.4])
    
    for tl in ax2.get_yticklabels():
        tl.set_color(Color2)
            
    if name in legendpars:             
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        legend = ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper center', fontsize=10, markerscale=0)
    
        frame = legend.get_frame()        
        frame.set_linewidth(0.4)   
        frame.set_edgecolor(ColorGrey)

    xmin = vals[0, 0] - 0.1 * (vals[-1, 0] - vals[0, 0])
    xmax = vals[-1, 0] + 0.1 * (vals[-1, 0] - vals[0, 0])
    ax2.set_xlim([xmin, xmax])

    [line.set_zorder(3) for line in ax1.lines]
    [line.set_zorder(3) for line in ax2.lines]
        
    fig.tight_layout()    
    fig.savefig(os.getcwd() + '//Robustness/puzzle_net_' + name + '.pdf')
    
    plt.close('all')     
    return
    
         
def find_vals(par):
           
    # a. load data
    filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
    f = h5py.File(filename, 'a')
    g = f.require_group('/sim')    
    sim = SimStruct()
         
    for x in ['pop_d','pop_n', 'pop_a', 'pop_P']:
                  
        exec('temp = g[\'' + x + '\']')         
        exec('sim.' + x + ' = np.copy(temp[...])')  
                   
    f.close()

    # b. calculations
    calc_sim(sim)
    t0 = 50
    
    P_all = sim.pop_P[t0, :]
    D_mean_all = np.mean(sim.pop_d[t0, :]*P_all)
    A_mean_all = np.mean(sim.pop_a[t0, :]*P_all)
    
    share_puzzle = sim.pop_puzzle[t0, :].sum()/par.N*100
    share_borrower = sim.pop_borrower[t0, :].sum()/par.N*100

    
    # c. return
    return np.array([share_puzzle, share_borrower, D_mean_all, A_mean_all])


def allfigs():
        
    path = os.getcwd() + '\\Robustness'
    if not os.path.exists(path):
        os.makedirs(path)
        
    files = glob.glob(path  + '//*.pdf')
    for f in files:
        os.remove(f)
        
    versions = [
    
                ['chi_lose', r'$\chi_{lose}$'],
    
                ['ra', r'$r_a$ (annual, fixed $r_d - r_a$)'],
                ['Gamma', r'$\Gamma$'],               

                ['sigma_xi', r'$\sigma_{\xi}$'],
                ['sigma_psi', r'$\sigma_{\psi}$'],
                ['u_ast', r'$u_{\ast}$'],
                ['pi_uu', r'$\pi_{u,u}$'],    

                ['lamb', r'$\lambda$'],                
                ['drd', r'$r_d$ (annual, fixed $r_a$)'],                                 

                ['eta', r'$\eta$'], 
                ['varphi', r'$\varphi$'],

                ['pi_x_lose_tot', r'$\pi^{lose}_x$'],
                ['pi_x_gain_tot', r'$\pi^{gain}_x$'],

                ['mu', r'$\mu$'],

                ['trans', r'$\varsigma$']   
                

                ]
                      
    legendpars = ['beta', 'chi_lose', 
                  'ra', 'trans', 'drd', 'pi_x_lose_tot']
    
    for version in versions:
        
        name = version[0] 
        name_latex = version[1]
        print(name)

        folders = glob.glob(os.getcwd() + '//' + name + '*')
        if name in ['ra'] and os.getcwd()+'\\random.h5' in folders:
            folders.remove(os.getcwd()+'\\random.h5')
        folders.append('baseline_med_0')
                       
        vals = np.zeros((len(folders), 5))    
        for i, folder in enumerate(folders):

            par = pickle.load(open( folder + '\\Data\\par.p', 'rb' ))
            
            if name in ['beta', 'Gamma']:            
                exec('vals[i, 0] = par.' + name + '**4')
            elif name in ['ra', 'drd']:            
                exec('vals[i, 0] = (1.0+par.' + name + ')**4-1')
            else:
                exec('vals[i, 0] = par.' + name)
                
            vals[i, 1:] = find_vals(par)
                        
        parbase = vals[-1, 0]   
        
        I = np.argsort(vals[:, 0])
        for j in xrange(len(vals[0, :])):
            vals[:, j] = vals[I, j]

        fig_robustness(parbase, name, name_latex,vals, legendpars)
     
    return

# main #
if __name__ == '__main__':
        
    allfigs()