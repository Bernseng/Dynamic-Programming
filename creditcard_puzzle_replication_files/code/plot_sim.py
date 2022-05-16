# -*- coding: utf-8 -*-
"""
Prints simulation figures and construct tables.

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

plt.rcParams.update({'font.size': 8,
                     'font.family': 'STIXGeneral',
                     'mathtext.fontset': 'stix'})

Color1 = (3.0/255.0,103.0/255.0,166.0/255.0)
Color2 = (242.0/255.0,62.0/255.0,46.0/255.0)
Color3 = (3.0/255.0,166.0/255.0,166.0/255.0)
Color4 = (242.0/255.0,131.0/255.0,68.0/255.0)
Color5 = (242.0/255.0,100.0/255.0,48.0/255.0)
ColorGrey = (65.0/255.0,68.0/255.0,81.0/255.0)

from tabulate import tabulate

class SimStruct:
    def __init__(self):
        pass


################
# CALCULATIONS #
################

def calc_sim(s):
    
    s.cutoff_base = 0.037
    s.cutoff_broad = 0.020
    s.cutoff_verybroad = 0.005    
    s.cutoff_narrow = 0.074
        
    cutoff = s.cutoff_base
    s.pop_puzzle = ((s.pop_a*s.pop_P > cutoff) & (s.pop_d*s.pop_P > cutoff))
    s.pop_borrower = (s.pop_puzzle == False) & (s.pop_d*s.pop_P-s.pop_a*s.pop_P > cutoff)
    s.pop_saver = (s.pop_puzzle == False) & (s.pop_a*s.pop_P-s.pop_d*s.pop_P > cutoff)
    s.pop_corner = (s.pop_puzzle == False) & (s.pop_borrower == False) & (s.pop_saver == False)

    cutoff = s.cutoff_broad
    s.pop_puzzle_broad = ((s.pop_a*s.pop_P > cutoff) & (s.pop_d*s.pop_P > cutoff))
    s.pop_borrower_broad = (s.pop_puzzle_broad == False) & (s.pop_d*s.pop_P-s.pop_a*s.pop_P > cutoff)
    s.pop_saver_broad = (s.pop_puzzle_broad == False) & (s.pop_a*s.pop_P-s.pop_d*s.pop_P > cutoff)
    s.pop_corner_broad = (s.pop_puzzle_broad == False) & (s.pop_borrower_broad == False) & (s.pop_saver_broad == False)
    
    cutoff = s.cutoff_verybroad
    s.pop_puzzle_verybroad = ((s.pop_a*s.pop_P > cutoff) & (s.pop_d*s.pop_P > cutoff))
    s.pop_borrower_verybroad = (s.pop_puzzle_verybroad == False) & (s.pop_d*s.pop_P-s.pop_a*s.pop_P > cutoff)
    s.pop_saver_verybroad = (s.pop_puzzle_verybroad == False) & (s.pop_a*s.pop_P-s.pop_d*s.pop_P > cutoff)
    s.pop_corner_verybroad = (s.pop_puzzle_verybroad == False) & (s.pop_borrower_verybroad == False) & (s.pop_saver_verybroad == False)    
    
    s.pop_puzzle_narrow = ((s.pop_a*s.pop_P > s.cutoff_narrow) & (s.pop_d*s.pop_P > s.cutoff_narrow))
    
    return   


########    
# TIME #
########

def fig_time_path_puzzle(p, s):
            
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)

    x_puzzle_mean = np.empty((p.simT-p.simBurnIn))
    x_puzzle_narrow_mean = np.empty((p.simT-p.simBurnIn))
    x_puzzle_broad_mean = np.empty((p.simT-p.simBurnIn)) 
    
    for t in xrange(p.simT-p.simBurnIn):
        
        x_puzzle_mean[t] = np.mean(np.array(s.pop_puzzle[t, :], dtype=float))*100
        x_puzzle_narrow_mean[t] = np.mean(np.array(s.pop_puzzle_narrow[t, :], dtype=float))*100
        x_puzzle_broad_mean[t] = np.mean(np.array(s.pop_puzzle_broad[t, :], dtype=float))*100         
            
    ax.plot(np.arange(p.simT-p.simBurnIn), x_puzzle_mean,
            lw=1, ls='-', marker='.', markersize=3,
            color=Color1,
            label=r'$a_t > ' + str(round(s.cutoff_base*100)/100) + '$' + ', $d_t > ' + str(round(s.cutoff_base*1000)/1000) + '$')
    
    ax.plot(np.arange(p.simT-p.simBurnIn), x_puzzle_narrow_mean,
            lw=1, ls='-', marker='.', markersize=3,
            color=Color2,
            label=r'$a_t > ' + str(round(s.cutoff_narrow*100)/100) + '$' + ', $d_t > ' + str(round(s.cutoff_narrow*1000)/1000) + '$')

    ax.plot(np.arange(p.simT-p.simBurnIn), x_puzzle_broad_mean,
            lw=1, ls='-', marker='.', markersize=3,
            color=Color3,
            label=r'$a_t > ' + str(round(s.cutoff_broad*100)/100) + '$' + ', $d_t > ' + str(round(s.cutoff_broad*1000)/1000) + '$')
            
    ax.set_xlabel(r'$t$', fontsize=10)
    ax.set_ylabel(r'percent', fontsize=10)
    ax.set_xlim([0, p.simT-p.simBurnIn-1])
    ax.set_ylim([np.fmax(np.min(x_puzzle_narrow_mean)-5, 0), np.fmin(np.max(x_puzzle_broad_mean)+5, 100)])

    legend = ax.legend(loc='best', fontsize=8)
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)

    [line.set_zorder(3) for line in ax.lines]
    fig.tight_layout()
    fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs//Sim//time_path_puzzle.pdf')
    
    plt.close('all') 
    return
        
        
def fig_time_path(p, s, x, name, latexname, mean, median, bounds):
        
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)

    if mean: 
        ax.plot(np.arange(p.simT-p.simBurnIn), np.mean(x, axis=1),
                lw=1.5, ls='-', marker='', markersize=3,
                color=Color1,
                label=r'mean')
                
    if median: 
        ax.plot(np.arange(p.simT-p.simBurnIn), np.median(x, axis=1),
                lw=1.5, ls='-', marker='', markersize=3,
                color=Color2,
                label=r'median')
            
    if bounds:
        ax.plot(np.arange(p.simT-p.simBurnIn), np.percentile(x, 25, axis=1),
                lw=1.5, ls='--', marker='', markersize=3,
                color=Color2)   
                
        ax.plot(np.arange(p.simT-p.simBurnIn), np.percentile(x, 75, axis=1),
                lw=1.5, ls='--', marker='', markersize=3,
                color=Color2)                   

    ax.set_xlabel(r'$t$', fontsize=10)
    ax.set_ylabel(r'$' + latexname + '_t$', fontsize=10)
    ax.set_xlim([0, p.simT-p.simBurnIn-1])
    
    if name in ['u', 'x']:
        ax.set_ylim([0, 20])
    
    legend = ax.legend(loc='best', fontsize=8)
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)

    [line.set_zorder(3) for line in ax.lines]
    fig.tight_layout()   
    fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/time_path_' + name + '.pdf')
    
    plt.close('all')
    return


########    
# DIST #
########

def fig_dist(p, sim, t, x, name, namelatex):

    # groups 
    x_all = x[t, :]
    x_puzzle = x[t, sim.pop_puzzle[t, :]]
    x_borrower = x[t, sim.pop_borrower[t, :]]
    x_saver = x[t, sim.pop_saver[t, :]]
    x_corner = x[t, sim.pop_corner[t, :]]
    
    if name in ['age']:
        bins_all = 401        
        bins = 401
    else:
        bins_all = 100
        bins = 50
        
    # all
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    plt.hist(x_all, bins_all, histtype="stepfilled", normed=True, color=Color1, alpha=.7, label='All');
    
    ax.set_xlabel(r'$' + namelatex + '$', fontsize=10)
   
    legend = ax.legend(loc='best', fontsize=8)
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)
        
    [line.set_zorder(3) for line in ax.lines]        
    fig.tight_layout()    
    fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/dist_all_' + name + '_' + str(t) + '.pdf')

    # puzzle vs. all
    if name not in ['age']:
            
        fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
        ax = fig.add_subplot(1,1,1)
            
        ax.spines["top"].set_linewidth(0.5) 
        ax.spines["bottom"].set_linewidth(0.5) 
        ax.spines["right"].set_linewidth(0.5)  
        ax.spines["left"].set_linewidth(0.5)  
        ax.tick_params(which='both', color=ColorGrey)
        ax.grid(True, zorder=0, color=ColorGrey)
        
        plt.hist(x_all, bins, histtype="stepfilled", normed=True, color=Color1, alpha=.7, label='All');
        plt.hist(x_puzzle, bins, histtype="stepfilled", normed=True, color=Color2, alpha=.7, label='Puzzle');
        
        ax.set_xlabel(r'$' + namelatex + '$', fontsize=10)
        legend = ax.legend(loc='best', fontsize=8)
        frame = legend.get_frame()        
        frame.set_linewidth(0.4)   
        frame.set_edgecolor(ColorGrey)

        [line.set_zorder(3) for line in ax.lines]        
        fig.tight_layout()
        fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/dist_puzzle_' + name + '_' + str(t) + '.pdf')
        
    # borrower vs. all
    if name not in ['age']:
        
        fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
        ax = fig.add_subplot(1,1,1)
            
        ax.spines["top"].set_linewidth(0.5) 
        ax.spines["bottom"].set_linewidth(0.5) 
        ax.spines["right"].set_linewidth(0.5)  
        ax.spines["left"].set_linewidth(0.5)  
        ax.tick_params(which='both', color=ColorGrey)
        ax.grid(True, zorder=0, color=ColorGrey)
        
        plt.hist(x_all, bins, histtype="stepfilled", normed=True, color=Color1, alpha=.7, label='All');
        plt.hist(x_borrower, bins, histtype="stepfilled", normed=True, color=Color2, alpha=.7, label='Borrower');
        
        ax.set_xlabel(r'$' + namelatex + '$', fontsize=10)
        legend = ax.legend(loc='best', fontsize=8)
        frame = legend.get_frame()        
        frame.set_linewidth(0.4)   
        frame.set_edgecolor(ColorGrey)
             
        [line.set_zorder(3) for line in ax.lines]        
        fig.tight_layout()        
        fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/dist_borrower_' + name + '_' + str(t) + '.pdf')
   
    # saver vs. all
    if name not in ['age']:    
        
        fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
        ax = fig.add_subplot(1,1,1)
            
        ax.spines["top"].set_linewidth(0.5) 
        ax.spines["bottom"].set_linewidth(0.5) 
        ax.spines["right"].set_linewidth(0.5)  
        ax.spines["left"].set_linewidth(0.5)  
        ax.tick_params(which='both', color=ColorGrey)
        ax.grid(True, zorder=0, color=ColorGrey)
        
        plt.hist(x_all, bins, histtype="stepfilled", normed=True, color=Color1, alpha=.7, label='All');
        plt.hist(x_saver, bins, histtype="stepfilled", normed=True, color=Color2, alpha=.7, label='Saver');
        
        ax.set_xlabel(r'$' + namelatex + '$', fontsize=10)
        legend = ax.legend(loc='best', fontsize=8)
        frame = legend.get_frame()        
        frame.set_linewidth(0.4)   
        frame.set_edgecolor(ColorGrey)
                
        [line.set_zorder(3) for line in ax.lines]        
        fig.tight_layout()        
        fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/dist_saver_' + name + '_' + str(t) + '.pdf')
        
    # corner vs. all
    if name not in ['age']:    
        
        fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
        ax = fig.add_subplot(1,1,1)
            
        ax.spines["top"].set_linewidth(0.5) 
        ax.spines["bottom"].set_linewidth(0.5) 
        ax.spines["right"].set_linewidth(0.5)  
        ax.spines["left"].set_linewidth(0.5)  
        ax.tick_params(which='both', color=ColorGrey)
        ax.grid(True, zorder=0, color=ColorGrey)
        
        plt.hist(x_all, bins, histtype="stepfilled", normed=True, color=Color1, alpha=.7, label='All');
        plt.hist(x_corner, bins, histtype="stepfilled", normed=True, color=Color2, alpha=.7, label='Corner');
        
        ax.set_xlabel(r'$' + namelatex + '$', fontsize=10)
        legend = ax.legend(loc='best', fontsize=8)
        frame = legend.get_frame()        
        frame.set_linewidth(0.4)   
        frame.set_edgecolor(ColorGrey)
        
        [line.set_zorder(3) for line in ax.lines]        
        fig.tight_layout()        
        fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/dist_corner_' + name + '_' + str(t) + '.pdf')
        
    plt.close('all')    
    return
    
    
########    
# AGE #
########

def fig_age(p, s, pop_age, x, name, namelatex):
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)

    x = np.ravel(x)
    x_mean = np.zeros(401)
    
    for age in xrange(401):
        I = np.ravel(pop_age) == age
        x_mean[age] = np.mean(x[I])
           
    ax.plot(np.arange(401)/4, x_mean,
            lw=1, ls='-', marker='.', markersize=3, color=Color1)
                
    ax.set_xlabel(r'age', fontsize=10)
    ax.set_ylabel(r'$' + namelatex + '$', fontsize=10)
    ax.set_xlim([0, p.simT-p.simBurnIn-1])

    [line.set_zorder(3) for line in ax.lines]        
    fig.tight_layout() 
    fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs//Sim//age_' + name + '.pdf')
    plt.close('all')
    return
   
def fig_age_num(p, s, pop_age):
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)

    x_num = np.zeros(401)
    for age in xrange(401):
        
        I = np.ravel(pop_age) == age
        x_num[age] = I.sum()
           
    ax.plot(np.arange(401)/4, x_num,
            lw=1, ls='-', marker='.', markersize=3, color=Color1)
                
    ax.set_xlabel(r'age', fontsize=10)
    ax.set_ylabel(r'#', fontsize=10)
    ax.set_xlim([0, p.simT-p.simBurnIn-1])

    [line.set_zorder(3) for line in ax.lines]        
    fig.tight_layout() 
    fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs//Sim//age_num.pdf')
    plt.close('all')
    return
    
    
################
# BEFORE/AFTER #
################
    
def fig_before_after_combined(p, sim, t0):
                    
    I = (sim.pop_puzzle[t0, :] == True) & (sim.pop_age[t0, :] > 0) & (sim.pop_puzzle[t0-1, :] == False) 
        
    k = 16
    mean = np.empty((2*k+1))
    median = np.empty((2*k+1))
    p25 = np.empty((2*k+1))
    p75 = np.empty((2*k+1))    
    
    names_latex_0 = [r'I. Puzzle share', 
                     r'II. Borrower share',
                     r'III. Saver share',
                     r'IV. Corner share']
                   
    names_latex_1 = [r'I. Liquid net worth - $n_t$', 
                     r'II. Unemployment - $u_t$', 
                     r'III. Permanent income - $P_t$', 
                     r'IV. Actual income - $Y_t$']

    for i in xrange(2):
        
        fig = plt.figure(frameon=False, figsize=(5.5, 4.5), dpi=800)   
        
        if i == 0:
            varsnow = ['puzzle', 'borrower', 'saver', 'corner']
            names_latex = names_latex_0
        else:
            varsnow = ['n', 'u', 'P', 'Y']
            names_latex = names_latex_1
            
        for j, name in enumerate(varsnow):
                    
            ax = fig.add_subplot(2,2,j+1) 
            
            ax.spines["top"].set_linewidth(0.5) 
            ax.spines["bottom"].set_linewidth(0.5) 
            ax.spines["right"].set_linewidth(0.5)  
            ax.spines["left"].set_linewidth(0.5)  
            ax.tick_params(which='both', color=ColorGrey)
            ax.grid(True, zorder=0, color=ColorGrey)
            
            exec("x = sim.pop_" + name)
                         
            for t in xrange(-k,k+1):
                
                mean[t+k] = np.mean(x[t0+t, I])
                
                if name not in ['P', 'Y', 'n']:   
                    mean[t+k] *= 100
                
                if name in ['P', 'Y', 'n']:
                
                    median[t+k] = np.median(x[t0+t, I])
                    p25[t+k] = np.percentile(x[t0+t, I], 25)
                    p75[t+k] = np.percentile(x[t0+t, I], 75)
            
            ax.plot(np.arange(-k,k+1), mean,
                    lw=1, ls='-', marker='o', markeredgecolor='none', markersize=2,
                    color=Color1,
                    label=r'mean')
            
            if name in ['P', 'Y', 'n']:               
                
                ax.plot(np.arange(-k,k+1), median,
                        lw=1, ls='-', marker='o', markeredgecolor='none', markersize=1.5,
                        color=Color2,
                        label=r'median')
                        
                ax.plot(np.arange(-k,k+1), p25,
                        lw=1, ls='--', marker='', markeredgecolor='none', markersize=1.5,
                        color=Color2,
                        label=r'25th percentile')
                        
                ax.plot(np.arange(-k,k+1), p75,
                        lw=1, ls='--', marker='', markeredgecolor='none', markersize=1.5,
                        color=Color2,
                        label=r'75th percentile')               
            else:
                
                ax.set_ylabel('percent', fontsize=7)        
    
            if name in ['P']:           
                legend = ax.legend(loc='best', fontsize=5)
                frame = legend.get_frame()        
                frame.set_linewidth(0.4)   
                frame.set_edgecolor(ColorGrey)
        
            ax.set_xlabel(r'$k$, quarters', fontsize=7)        
            ax.grid(True)
            ax.set_title(names_latex[j], fontsize=9)
            ax.tick_params(axis='both', which='major', labelsize=6)
            
            plt.axvline(x=0, ymin=0, ymax=100, linewidth=1, color='black', alpha=0.5)
            
        [line.set_zorder(3) for line in ax.lines]        
        fig.tight_layout()       
        if i == 0:
            fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/time_before_after_combined_groups.pdf')
        else:
            fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/time_before_after_combined_variables.pdf')

    plt.close('all') 
    return
            
            
def fig_before_after(p, sim, t0, x, name, justmean):
    
    fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)

    I = (sim.pop_puzzle[t0, :] == True) & (sim.pop_age[t0, :] > 0) & (sim.pop_puzzle[t0-1, :] == False) 
        
    k = 16
    mean = np.empty((2*k+1))
    median = np.empty((2*k+1))
    p25 = np.empty((2*k+1))
    p75 = np.empty((2*k+1))    
    
    for t in xrange(-k,k+1):
        
        mean[t+k] = np.mean(x[t0+t, I])
        median[t+k] = np.median(x[t0+t, I])
        p25[t+k] = np.percentile(x[t0+t, I], 25)
        p75[t+k] = np.percentile(x[t0+t, I], 75)
    
    ax.plot(np.arange(-k,k+1), mean,
            lw=1.5, ls='-', marker='o', markeredgecolor='none', markersize=3,
            color=Color1,
            label=r'mean')

    if justmean == 0:
        
        ax.plot(np.arange(-k,k+1), median,
                lw=1.5, ls='-', marker='', markeredgecolor='none', markersize=2,
                color=Color2,
                label=r'median')
                
        ax.plot(np.arange(-k,k+1), p25,
                lw=1.5, ls='--', marker='', markeredgecolor='none',markersize=2,
                color=Color2)
                
        ax.plot(np.arange(-k,k+1), p75,
                lw=1.5, ls='--', marker='', markeredgecolor='none', markersize=2,
                color=Color2)            

    ax.set_xlabel(r'$k$, quarters', fontsize=10)
    ax.grid(True)
    
    if justmean == 0:            
        legend = ax.legend(loc='best', fontsize=8)
        frame = legend.get_frame()        
        frame.set_linewidth(0.4)   
        frame.set_edgecolor(ColorGrey)
    
    plt.axvline(x=0, ymin=0, ymax=100, linewidth=1, color='black', alpha=0.5)
            
    [line.set_zorder(3) for line in ax.lines]
    fig.tight_layout()    
    fig.savefig(os.getcwd() + '\\' + p.name + '//Graphs/Sim/time_before_after_' + name + '_' + str(t0) + '.pdf')
    
    plt.close('all')    
    return    
    
    
#########    
# TABLE #
#########

def create_standard_block(table, sim, x, data, t0, name, do_all = 1):
    
    P_puzzle = sim.pop_P[t0, sim.pop_puzzle[t0, :]]
    if sim.pop_puzzle[t0, :].sum() > 0:
        x_mean_puzzle = np.mean(x[t0, sim.pop_puzzle[t0, :]]*P_puzzle)
    else:
        x_mean_puzzle = ''
        
    P_borrower = sim.pop_P[t0, sim.pop_borrower[t0, :]]
    if sim.pop_borrower[t0, :].sum() > 0:     
        x_mean_borrower = np.mean(x[t0, sim.pop_borrower[t0, :]]*P_borrower)
    else:
        x_mean_borrower = ''
        
    P_saver = sim.pop_P[t0, sim.pop_saver[t0, :]]
    if sim.pop_saver[t0, :].sum() > 0:      
        x_mean_saver = np.mean(x[t0, sim.pop_saver[t0, :]]*P_saver)
    else:
        x_mean_saver = ''
        
    P_corner = sim.pop_P[t0, sim.pop_corner[t0, :]]
    if sim.pop_corner[t0, :].sum() > 0:    
        x_mean_corner = np.mean(x[t0, sim.pop_corner[t0, :]]*P_corner)
    else:
        x_mean_corner = ''
        
    P_all = sim.pop_P[t0, :]
    x_mean_all = np.mean(x[t0, :]*P_all)
         
    row = [name, x_mean_puzzle, 
                 x_mean_borrower,
                 x_mean_saver,
                 x_mean_corner,
                 x_mean_all, data[0]]
    
    table.append(row)
    
    if do_all == 1:
        percs = [5, 15, 25, 50, 75, 85, 95]
    else:
        percs = [25, 50, 75]

    for i, percnow in enumerate(percs):
        
        P_puzzle = sim.pop_P[t0, sim.pop_puzzle[t0, :]]
        if sim.pop_puzzle[t0, :].sum() > 0:
            x_perc_puzzle = np.percentile(x[t0, sim.pop_puzzle[t0, :]]*P_puzzle,percnow)
        else:
            x_perc_puzzle = ''
            
        P_borrower = sim.pop_P[t0, sim.pop_borrower[t0, :]]
        if sim.pop_borrower[t0, :].sum() > 0:   
            x_perc_borrower = np.percentile(x[t0, sim.pop_borrower[t0, :]]*P_borrower,percnow)
        else:
            x_perc_borrower = ''
            
        P_saver = sim.pop_P[t0, sim.pop_saver[t0, :]]
        if sim.pop_saver[t0, :].sum() > 0:   
            x_perc_saver = np.percentile(x[t0, sim.pop_saver[t0, :]]*P_saver,percnow)
        else:
            x_perc_saver = ''
            
        P_corner = sim.pop_P[t0, sim.pop_corner[t0, :]]        
        if sim.pop_corner[t0, :].sum() > 0:   
            x_perc_corner = np.percentile(x[t0, sim.pop_corner[t0, :]]*P_corner,percnow)
        else:
            x_perc_corner = ''
            
        P_all = sim.pop_P[t0, :]
        x_perc_all = np.percentile(x[t0, :]*P_all,percnow)
             
        row = ['\\hspace{3mm}' + str(percnow) + 'th percentile', 
                   x_perc_puzzle, 
                   x_perc_borrower,
                   x_perc_saver,
                   x_perc_corner,
                   x_perc_all, data[1+i]]
    
        table.append(row)
    
    return
    
def print_table(par, sim, t0):
        
    table = []

    # a. share row    
    simshape = sim.pop_puzzle.shape
    N = simshape[1];
    row = ['Share',sim.pop_puzzle[t0, :].sum()/N*100, 
                   sim.pop_borrower[t0, :].sum()/N*100, 
                   sim.pop_saver[t0, :].sum()/N*100,
                   sim.pop_corner[t0, :].sum()/N*100,
                   100]                                  
    table.append(row)  
    
    # b. other rows
    create_standard_block(table, sim, sim.pop_d, par.d_data, t0, 'Credit Card Debt, $D_t$ (mean)')
    create_standard_block(table, sim, sim.pop_a, par.a_data, t0, 'Liquid assets, $A_t$ (mean)')
    create_standard_block(table, sim, sim.pop_n, par.n_data, t0, 'Liquid net worth, $N_t$ (mean)')
    
    # c. print table in plain text
    header = ['Puzzle', 'Borrower', 'Saver', 'Corner', 'All'] 
    with open(os.getcwd() + '\\' + par.name + "\\table.txt", "w") as text_file:
        text_file.write(par.tableline)
        text_file.write("\n")
        text_file.write(tabulate(table, headers=header, 
                        tablefmt="pipe",
                        floatfmt=".2f"))    
        text_file.write("\n\n")

    # c. print table in latex    
    with open(os.getcwd() + '\\' + par.name + "//table_latex.txt", "w") as text_file:
       
        headrow = ''
        headrow += " & \\multicolumn{1}{c}{Puzzle}"
        headrow += " & \\multicolumn{1}{c}{Borrower}"
        headrow += " & \\multicolumn{1}{c}{Saver}"
        headrow += " & \\multicolumn{1}{c}{Corner}"
        headrow += " & \\multicolumn{1}{c}{All}"    
        headrow += " & \\multicolumn{1}{c}{Data}"  
        headrow += " \\" + "\\" + "\n"         
        text_file.write(headrow) 
        text_file.write("\\midrule \n") 
        text_file.write("\\addlinespace \n")
        text_file.write(" & \multicolumn{6}{c}{\\textit{percent}}")
        text_file.write(" \\" + "\\ \n")
        text_file.write("\\addlinespace \n")

        for i, row in enumerate(table):   
                
            string = row[0]
            if i in [0]:
                for column in row[1:]:
                    string += ' & '
                    if column != '':
                        string += '{:3.1f} '.format(column) 
            else:
                for j, column in enumerate(row[1:]):
                    string += ' & '
                    if column != '':
                        string += '{:3.2f} '.format(column)
                        
            string += ' \\' + '\\' + '\n'
            text_file.write(string) 
            
            if i == 0: 
                text_file.write("\\addlinespace \n")  
                text_file.write("\\midrule \n")   
                text_file.write("\\addlinespace \n")
                text_file.write(" & \multicolumn{6}{c}{\\textit{relative to mean quarterly income}}")
                text_file.write(" \\" + "\\ \n")
                text_file.write("\\addlinespace \n")
            elif i in [1,8,9,16,17]: 
                text_file.write("\\addlinespace \n\n")
                
    return
    
def print_robustness_table():
    
    path = os.getcwd() + '\\Robustness'
    if not os.path.exists(path):
        os.makedirs(path)
        
    t0 = 50
    names_tuple = [
    
                   ['baseline_0','Baseline'],

                   ['baseline_chi_lose_0','$\\chi_{lose} = 1.0$'], 
                   ['baseline_chi_lose_1','$\\chi_{lose} = 2.0$'],
                   ['baseline_chi_lose_2','$\\chi_{lose} = 6.0$'],
                    
                   ['baseline_pi_x_lose_tot_0','$\\pi_{\\ast}^{lose} = 0.005$'], 
                   ['baseline_pi_x_lose_tot_1','$\\pi_{\\ast}^{lose} = 0.010$'],
                   ['baseline_pi_x_lose_tot_2','$\\pi_{\\ast}^{lose} = 0.040$'],

                   ['baseline_pi_x_gain_tot_0','$\\pi_{\\ast}^{gain} = 0.03$'], 
                   ['baseline_pi_x_gain_tot_1','$\\pi_{\\ast}^{gain} = 0.12$'],
                   ['baseline_pi_x_gain_tot_2','$\\pi_{\\ast}^{gain} = 0.24$'],
                    
                   ['baseline_u_ast_0','$u_{\\ast} = 0.04$'], 
                   ['baseline_u_ast_1','$u_{\\ast} = 0.05$'],
                   ['baseline_u_ast_2','$u_{\\ast} = 0.06$']
                    

                  ]
                  
    with open(os.getcwd() + '\\Robustness\\' + ' table_robustness.txt', 'w') as text_file:
       
        headrow = ''
        headrow += " & \\multicolumn{1}{c}{Puzzle}"
        headrow += " & \\multicolumn{1}{c}{Borrower}"
        headrow += " & \\multicolumn{1}{c}{Saver}"
        headrow += " & \\multicolumn{1}{c}{Corner}" 
        headrow += " & \\multicolumn{1}{c}{$A_t$}"
        headrow += " & \\multicolumn{1}{c}{$D_t$}"
        headrow += " & \\multicolumn{1}{c}{$N_t$}"
        headrow += " \\" + "\\" + "\n"         
        text_file.write(headrow) 
        text_file.write("\\midrule \n") 
        text_file.write("\\addlinespace \n")
        text_file.write(" & \multicolumn{4}{c}{\\textit{percent}} & \multicolumn{3}{c}{\\textit{mean}}")
        text_file.write(" \\" + "\\ \n")
        text_file.write("\\addlinespace \n")

        sim = SimStruct()
        for name_tuple in names_tuple:
            
            name = name_tuple[0]
            name_latex = name_tuple[1] 
            
            # i. load variables            
            par = pickle.load(open(os.getcwd() + '\\' + name + '\\Data\\par.p', 'rb' ))                    
            
            for x in ['pop_d', 'pop_n', 'pop_a', 'pop_P']:
                      
                ini = 1              
                for i in xrange(0,par.nbeta):
                    for j in xrange(0,par.nrho):
            
                        par.het = '' + str(i) + '_' + str(j)
                    
                        filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
                        f = h5py.File(filename, 'a')
                        g = f.require_group('/sim')    
                        
                        exec('temp = g[\'' + x + '\']')                  
                        if ini == 1:             
                            exec('sim.' + x + ' = np.copy(temp[...])')  
                            ini = 0
                        else:
                            exec('sim.' + x + ' = np.concatenate((sim.' + x + ',np.copy(temp[...])),axis=1)') 
                        exec('sim.' + x + '.shape')
                                 
                        f.close()
                        
            calc_sim(sim)
            
            # ii. calculate shares
            simshape = sim.pop_puzzle.shape
            N = simshape[1]
            row = []
            row.append(sim.pop_puzzle[t0, :].sum()/N*100)
            row.append(sim.pop_borrower[t0, :].sum()/N*100)
            row.append(sim.pop_saver[t0, :].sum()/N*100)
            row.append(sim.pop_corner[t0, :].sum()/N*100)
            
            # iii. calculate mean
            P_all = sim.pop_P[t0, :]
            row.append(np.mean(sim.pop_d[t0, :]*P_all))
            row.append(np.mean(sim.pop_a[t0, :]*P_all))
            row.append(np.mean(sim.pop_n[t0, :]*P_all))
                        
            # iv. print row
            string = name_latex
            for i, column in enumerate(row):
                if i < 4:
                    string += ' & {:3.1f} '.format(column)
                else:
                    string += ' & {:3.2f} '.format(column)
            string += ' \\' + '\\' + '\n'
            text_file.write(string) 
            
            if name in ['baseline_0','baseline_chi_lose_2','baseline_pi_x_lose_tot_2','baseline_pi_x_gain_tot_2','baseline_u_ast_2']:
                text_file.write('\\addlinespace\n\n')
                    
            # iv. alternative cutoffs
            if name in ['baseline_0']:
                
                # broad
                row = []
                row.append(sim.pop_puzzle_broad[t0, :].sum()/N*100)
                row.append(sim.pop_borrower_broad[t0, :].sum()/N*100)
                row.append(sim.pop_saver_broad[t0, :].sum()/N*100)
                row.append(sim.pop_corner_broad[t0, :].sum()/N*100)
                
                # iii. calculate mean
                P_all = sim.pop_P[t0, :]
                row.append(np.mean(sim.pop_d[t0, :]*P_all))
                row.append(np.mean(sim.pop_a[t0, :]*P_all))
                row.append(np.mean(sim.pop_n[t0, :]*P_all))
                
                string = '\hspace{1mm}cut-off = 0.025'
                for i, column in enumerate(row):
                    if i < 4:
                        string += ' & {:3.1f} '.format(column)
                    else:
                        string += ' & {:3.2f} '.format(column)
                string += ' \\' + '\\' + '\n'
                text_file.write(string)             
            
                # very broad
                row = []
                row.append(sim.pop_puzzle_verybroad[t0, :].sum()/N*100)
                row.append(sim.pop_borrower_verybroad[t0, :].sum()/N*100)
                row.append(sim.pop_saver_verybroad[t0, :].sum()/N*100)
                row.append(sim.pop_corner_verybroad[t0, :].sum()/N*100)
                
                # iii. calculate mean
                P_all = sim.pop_P[t0, :]
                row.append(np.mean(sim.pop_d[t0, :]*P_all))
                row.append(np.mean(sim.pop_a[t0, :]*P_all))
                row.append(np.mean(sim.pop_n[t0, :]*P_all))
                
                string = '\hspace{1mm}cut-off = 0.005'
                for i, column in enumerate(row):
                    if i < 4:
                        string += ' & {:3.1f} '.format(column)
                    else:
                        string += ' & {:3.2f} '.format(column)
                string += ' \\' + '\\' + '\n'
                text_file.write(string) 
                
                text_file.write('\\addlinespace\n')
                
    return
        
###########    
# MOMENTS #
###########
        
def moments(par,sim,t0):
    
    simshape = sim.pop_puzzle.shape
    N = simshape[1];
    par.sim_puzzle_share = sim.pop_puzzle[t0, :].sum()/N*100
    par.sim_borrower_share = sim.pop_borrower[t0, :].sum()/N*100
    par.sim_saver_share = sim.pop_saver[t0, :].sum()/N*100
    par.sim_corner_share = sim.pop_corner[t0, :].sum()/N*100
    
    par.table = []
    create_standard_block(par.table, sim, sim.pop_d, par.d_data, t0, '')
    create_standard_block(par.table, sim, sim.pop_a, par.a_data, t0, '')
    create_standard_block(par.table, sim, sim.pop_n, par.n_data, t0, '')
      
    base = 1
    d_p5  = par.table[base][5]
    d_p15 = par.table[base+1][5]
    d_p25 = par.table[base+2][5]
    d_p50 = par.table[base+3][5]
    d_p75 = par.table[base+4][5]
    d_p85 = par.table[base+5][5]
    
    base = 9
    a_p5  = par.table[base][5]
    a_p15 = par.table[base+1][5]
    a_p25 = par.table[base+2][5]
    a_p50 = par.table[base+3][5]
    a_p75 = par.table[base+4][5]
    a_p85 = par.table[base+5][5]
    
    base = 17
    n_p5  = par.table[base][5]
    n_p15 = par.table[base+1][5]
    n_p25 = par.table[base+2][5]
    n_p50 = par.table[base+3][5]
    n_p75 = par.table[base+4][5]
    n_p85 = par.table[base+5][5]
    
    moms = [
            [d_p5, par.d_data[1]],
            [d_p15, par.d_data[2]],
            [d_p25, par.d_data[3]],
            [d_p50, par.d_data[4]],
            [d_p75, par.d_data[5]],
            [d_p85, par.d_data[6]],
            [a_p5,  par.a_data[1]],
            [a_p15, par.a_data[2]],
            [a_p25, par.a_data[3]],
            [a_p50, par.a_data[4]],
            [a_p75, par.a_data[5]],
            [a_p85, par.a_data[6]],
           ]
           
    print('')
    val = 0
    for i, mom in enumerate(moms):
        cont = ((mom[1] - mom[0])/100)**2.0;
        val += cont
        print("mom {:2d}: {:5.1f} vs. {:5.1f} [{:3.4f}]".format(i,100*mom[0],100*mom[1],100*100*cont))
    print('')
    
    return val

    
##########    
# SHARES #
##########
    
def shares(par,t0):
    
    p = np.zeros((par.nbeta,par.nrho))
    b = np.zeros((par.nbeta,par.nrho))
    s = np.zeros((par.nbeta,par.nrho))
    c = np.zeros((par.nbeta,par.nrho))
    for i in xrange(0,par.nbeta):
        for j in xrange(0,par.nrho):
                
            par.het = '' + str(i) + '_' + str(j)
            sim = SimStruct()
                
            for x in ['pop_d', 'pop_a', 'pop_P']:
              
                filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
                f = h5py.File(filename, 'a')
                g = f.require_group('/sim')    
                
                exec('temp = g[\'' + x + '\']')                        
                exec('sim.' + x + ' = np.copy(temp[...])')  
                         
                f.close()
        
            calc_sim(sim)

            simshape = sim.pop_puzzle.shape
            N = simshape[1];
            p[i,j] = sim.pop_puzzle[t0, :].sum()/N*100
            b[i,j] = sim.pop_borrower[t0, :].sum()/N*100
            s[i,j] = sim.pop_saver[t0, :].sum()/N*100
            c[i,j] = sim.pop_corner[t0, :].sum()/N*100
            
    for x, name in zip([p, b, s, c],['p','b','s','c']):
        
        fig = plt.figure(frameon=False, figsize=(5.5, 3.9), dpi=800)   
        ax = fig.add_subplot(1,1,1)
            
        ax.spines["top"].set_linewidth(0.5) 
        ax.spines["bottom"].set_linewidth(0.5) 
        ax.spines["right"].set_linewidth(0.5)  
        ax.spines["left"].set_linewidth(0.5)  
        ax.tick_params(which='both', color=ColorGrey)
        ax.grid(True, zorder=0, color=ColorGrey)
            
        ax.plot(par.rho_mat[0,:], x[0,:],
                lw=1.5, ls='-', marker='.', markersize=3,
                color=Color1,
                label=r'$\beta = {:5.3f}$'.format(par.beta_mat[0,0]))
        
        if par.nbeta > 2:
            ax.plot(par.rho_mat[2,:], x[2,:],
                    lw=2, ls='--', marker='.', markersize=5,
                    color=Color2,
                    label=r'$\beta = {:5.3f}$'.format(par.beta_mat[2,0]))

        if par.nbeta > 4:
            ax.plot(par.rho_mat[4,:], x[4,:],
                    lw=2, ls=':', marker='.', markersize=5,
                    color=Color3,
                    label=r'$\beta = {:5.3f}$'.format(par.beta_mat[4,0]))
                
        ax.set_xlabel(r'$\rho$', fontsize=11)
        ax.set_ylabel(r'percent', fontsize=11)
        
        if name in ['p']:
            legend = ax.legend(loc='best', fontsize=10)
            frame = legend.get_frame()        
            frame.set_linewidth(0.4)   
            frame.set_edgecolor(ColorGrey)
         
        [line.set_zorder(3) for line in ax.lines]
        fig.tight_layout()
        fig.savefig(os.getcwd() + '\\' + par.name + '\\Graphs\\Sim\\shares_prefs_' + name + '.pdf')
    
    plt.close('all')
    return


########    
# MENU #
########
          
def allfigs(par, early_quit = True):

    print("plot_sim.py")
    files = glob.glob(os.getcwd() + '\\' + par.name + '\\Graphs\\Sim\\*.pdf')
    for f in files:
        os.remove(f)
        
    # set cross-section time
    t0 = 50
       
    # a. shares for each preference group
    if par.nbeta*par.nrho > 1:
        shares(par,t0)
    
    # b. load all data  
    t1 = time.time()
    sim = SimStruct()
         
    for x in ['pop_age', 'pop_db', 'pop_u', 'pop_x', 'pop_db', 'pop_nb', 
              'pop_d', 'pop_c', 'pop_n', 'pop_a', 
              'pop_death', 'pop_psi', 'pop_xi',
              'pop_P', 'pop_Y',]:
              
        ini = 1              
        for i in xrange(0,par.nbeta):
            for j in xrange(0,par.nrho):
            
                par.het = '' + str(i) + '_' + str(j)
            
                filename = os.getcwd() + '\\' + par.name + '\\Data\\Data_' + str(par.het) + '.h5'
                f = h5py.File(filename, 'a')
                g = f.require_group('/sim')    
                
                exec('temp = g[\'' + x + '\']')                  
                if ini == 1:             
                    exec('sim.' + x + ' = np.copy(temp[...])')  
                    ini = 0
                else:
                    exec('sim.' + x + ' = np.concatenate((sim.' + x + ',np.copy(temp[...])),axis=1)') 
                exec('sim.' + x + '.shape')
                         
                f.close()

    print(" - data loaded ({:3.1f})".format(time.time()-t1))
    
    # c. calculations
    t1 = time.time()
    calc_sim(sim)
    print(" - calculations done ({:3.1f})".format(time.time()-t1))

    # d. table
    t1 = time.time()
    print_table(par, sim, t0)
    print(" - table printed ({:3.1f})".format(time.time()-t1))

    # e. moments
    moments(par, sim, t0)
        
    if early_quit:
        return
    
    # f. print figures
    t1 = time.time()        
    fig_before_after_combined(par, sim, t0)
    fig_before_after(par, sim, t0, sim.pop_puzzle*100, 'puzzle', 1)
    fig_before_after(par, sim, t0, sim.pop_borrower*100, 'borrower', 1)
    fig_before_after(par, sim, t0, sim.pop_saver*100, 'saver', 1)
    fig_before_after(par, sim, t0, sim.pop_corner*100, 'corner', 1)
    fig_before_after(par, sim, t0, sim.pop_a, 'a', 0)
    fig_before_after(par, sim, t0, sim.pop_d, 'd', 0)    
    fig_before_after(par, sim, t0, sim.pop_n, 'n', 0)  
    fig_before_after(par, sim, t0, sim.pop_c, 'c', 0)
    fig_before_after(par, sim, t0, sim.pop_u*100, 'u', 1)
    fig_before_after(par, sim, t0, sim.pop_x*100, 'x', 1)
    fig_before_after(par, sim, t0, sim.pop_psi, 'psi', 0)  
    fig_before_after(par, sim, t0, sim.pop_xi, 'xi', 0)
    fig_before_after(par, sim, t0, sim.pop_P, 'P', 0)
    fig_before_after(par, sim, t0, sim.pop_Y, 'Y', 0)  
        
    print(" - figs: before/after ({:3.1f})".format(time.time()-t1))
    
    t1 = time.time()
    fig_time_path_puzzle(par, sim) 
    
    fig_time_path(par, sim, sim.pop_n, 'n', 'n', 1, 1, 1) 
    fig_time_path(par, sim, sim.pop_a, 'a', 'a', 1, 1, 1)
    fig_time_path(par, sim, sim.pop_d, 'd', 'd', 1, 1, 1)     
    fig_time_path(par, sim, sim.pop_c, 'c', 'c', 1, 1, 1)
    
    fig_time_path(par, sim, sim.pop_u*100, 'u', 'u', 1, 0, 0)
    fig_time_path(par, sim, sim.pop_x*100, 'x', 'x', 1, 0, 0)
    fig_time_path(par, sim, sim.pop_age/4, 'age', 'age', 1, 1, 1)
    fig_time_path(par, sim, sim.pop_psi, 'psi', '\psi', 1, 1, 1)
    fig_time_path(par, sim, sim.pop_xi, 'xi', '\\xi', 1, 1, 1)

    fig_time_path(par, sim, sim.pop_P, 'P', 'P', 1, 1, 1)     
    fig_time_path(par, sim, sim.pop_Y, 'Y', 'Y', 1, 1, 1)
    
    print(" - figs: time path ({:3.1f})".format(time.time()-t1))
    
    t1 = time.time()
    fig_dist(par, sim, t0, sim.pop_age/4, 'age', 'age')    
    fig_dist(par, sim, t0, sim.pop_n, 'n', 'n_t')
    fig_dist(par, sim, t0, sim.pop_d, 'd', 'd_t')
    fig_dist(par, sim, t0, sim.pop_a, 'a', 'a_t')
    fig_dist(par, sim, t0, sim.pop_P, 'P', 'P_t')
    fig_dist(par, sim, t0, sim.pop_Y, 'Y', 'Y_t')
    
    print(" - figs: distributions ({:3.1f})".format(time.time()-t1))
                
    t1 = time.time()
    fig_age_num(par, sim, sim.pop_age)    
    fig_age(par, sim, sim.pop_age, sim.pop_puzzle, 'puzzle', 'puzzle')
    fig_age(par, sim, sim.pop_age, sim.pop_P, 'P', 'P_t')
    
    print(" - figs: age ({:3.1f})".format(time.time()-t1))
    return

        
# main #
if __name__ == '__main__':
    
    cases = [] #['baseline_0']
    print_robustness_table()

    for case in cases:
        
        # a. load data
        par = pickle.load(open( case + '\\Data\\par.p', 'rb' ))
        
        # b. print figures
        allfigs(par, False)