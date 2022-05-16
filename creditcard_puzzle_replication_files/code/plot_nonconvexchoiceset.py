# -*- coding: utf-8 -*-
"""
Plots the choice set.

Version: 1.0.
@author: Jeppe Druedahl, 2017.
"""

from __future__ import division, print_function, absolute_import

import numpy as np
import os
import glob

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
    
def fig_choiceset(db, nb, eta, varphi, name, leq):

    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
    
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    # 1. dt and c bounds
    d_min = np.fmax(-nb, 0)
    d_max = np.fmax(db, eta*nb+varphi)

    def c_max_d(db, nb, d, eta, varphi):
        if d <= db:
            return nb + d
        else:
            if eta > 0:
                return nb + np.fmin(d, 1/eta*(varphi-d))
            else:
                return nb + d
            
    # 2. c bounds given d
    d_c_max = np.fmax(db, varphi/(1+eta))
    c_d_min = c_max_d(db, nb, d_min, eta, varphi)
    c_d_max = c_max_d(db, nb, d_max, eta, varphi)
    c_d_db_l = c_max_d(db, nb, d_c_max, eta, varphi)
    c_d_db_r = c_max_d(db, nb, d_c_max+10**(-8), eta, varphi)
    c_max = nb + np.fmax(db, varphi/(1+eta))
    
    if d_max > db:
        assert np.allclose(c_d_max, 0), (d_max, db, c_d_max)

    # 3. line segments
    ax.plot([d_min, d_min], [0.0, c_d_min], '-o', color='black', markersize=4, alpha=0.50)
    ax.plot([d_min, d_c_max], np.zeros(2), '-o', color='black', markersize=4, alpha=0.50)
    ax.plot([d_min, d_c_max], [c_d_min, c_d_db_l], '-o', color='black', markersize=4, alpha=0.50)
    if d_max > d_c_max:
        ax.plot([d_c_max, d_max], np.zeros(2), '-o', color='black', markersize=4, alpha=0.50)
        ax.plot([d_c_max, d_max], [c_d_db_r, 0], '-o', color='black', markersize=4, alpha=0.50)
        if np.allclose(c_d_db_r, c_d_db_l) == 0:
            ax.plot([d_c_max, d_c_max], [c_d_db_r, c_d_db_l], '-o', color='black', markersize=4, alpha=0.50)
  
    else:
        ax.plot([d_c_max, d_c_max], [0, c_d_db_l], '-o', color='black', markersize=4, alpha=0.50)
        
    # 4. annotations
    if db == 0:
        ax.text(db, -0.1*c_max, r'$\bar{d}_t = d_{min}$', fontsize=10, color='black')
        ax.text(d_max, -0.1*c_max, r'$d_{max}$', fontsize=10, color='black')        
    elif d_c_max < d_max:
        ax.text(d_min, -0.1*c_max, r'$d_{min}$', fontsize=10, color='black')
        ax.text(d_max, -0.1*c_max, r'$d_{max}$', fontsize=10, color='black')
        ax.text(db, -0.1*c_max, r'$\bar{d}_t$', fontsize=10, color='black')
    else:
        ax.text(d_min, -0.1*c_max, r'$d_{min}$', fontsize=10, color='black')        
        ax.text(db, -0.1*c_max, r'$\bar{d}_t = d_{max}$', fontsize=10, color='black')
        
    if leq:
        ax.text(d_min+0.05, 0.1, r"$c_t \leq \,\ \bar{n}_t + d_t$", 
                fontsize=10, color="black")    
        ax.text(db+0.05, 0.1, 
                r'$c_t \leq \,\ \bar{n}_{t}+\min\left\{ d_{t}' 
                r',\,\frac{1}{\eta}\left(\varphi-d_{t}\right)\right\} $', 
                fontsize=10, color="black")

    # 5. title, labels and grid
    ax.set_title(r'$\bar{d}_t = ' + str(db) + r'$, ' \
                 + r'$\bar{n}_t = ' + str(nb) + r'$, ' \
                 + r'$\eta = ' + str(eta) + r'$, ' \
                 + r'$\varphi = ' + str(varphi) + r'$', fontsize=10)        
    
    ax.set_ylabel(r'$c_t$', fontsize=10)
    ax.set_xlabel(r'$d_t$', fontsize=10)

    ax.set_ylim(-0.2*c_max, c_max*1.2)
    ax.set_xlim(-0.2*d_max+d_min, d_max*1.2)
        
    # 6. fills
    if d_c_max > 0:
        ax.fill_between([d_min, d_c_max], [0, 0], [0, c_d_min], color=Color1)
        ax.fill_between([d_min, d_c_max], [0, 0], [c_d_min, c_d_db_l], color=Color1)

    if d_max > d_c_max:
        ax.fill_between([d_c_max, d_max], [0, 0], [c_d_db_r, 0], color=Color2)

    # 7. save   
    [line.set_zorder(3) for line in ax.lines]
    fig.tight_layout()
    fig.savefig('ChoiceSetGraphs\\choiceset' + name + '.pdf')
    plt.close('all') 
    return
    
if __name__ == '__main__':

    path = os.getcwd() + '\\ChoiceSetGraphs'
    if not os.path.exists(path):
        os.makedirs(path)
     
    files = glob.glob(path + '//*.pdf')
    for f in files:
        os.remove(f)
           
    eta = 1.0
    varphi = 0.5
    nb = 1.0
    
    db = 0.0
    fig_choiceset(db, nb, eta, varphi, '_zero_db', False)

    db = varphi/(1+eta)
    fig_choiceset(db, nb, eta, varphi, '_breakeven_db', False)
    
    db = 0.5
    fig_choiceset(db, nb, eta, varphi, '_low_db', True)
    
    db = 0.5
    fig_choiceset(db, nb, 0.0, varphi, '_low_db_zero_eta', False)
    
    db = 2.0
    fig_choiceset(db, nb, eta, varphi, '_high_db', False)

    db = 0.0
    nb = 0.5
    fig_choiceset(db, nb, eta, varphi, '_zero_db_low_nb', False)
    
    db = 0.5
    nb = -0.1
    fig_choiceset(db, nb, eta, varphi, '_low_db_neg_nb', False)