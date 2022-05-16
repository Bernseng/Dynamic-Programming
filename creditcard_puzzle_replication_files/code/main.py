# -*- coding: utf-8 -*-
"""
*MAIN SCRIPT*

Version: 10.
@author: Jeppe Druedahl, 2017.
"""

# imports #
from __future__ import division, print_function, absolute_import
import os
clear = lambda: os.system('cls')
clear()

import time
import glob
import h5py
from scipy.optimize import minimize

import numpy as np
import cPickle as pickle

import settings # used in an exec()
import grid
import solve
import sim
import plot_sol
import plot_sim
import plot_all
import plot_robustness
import plot_welfare


############################
# VALUE FUNCTION ITERATION #
############################

def vfi_run(par):
            
    counter = 0
    for i in xrange(0,par.nbeta):
        for j in xrange(0,par.nrho):
            
            par.het = '' + str(i) + '_' + str(j)
            par.beta = par.beta_mat[i,j]
            par.rho = par.rho_mat[i,j]  

            if par.nbeta*par.nrho > 1:
                counter += 1
                print('\n#het ' + str(counter) + ' of ' + str(par.nbeta*par.nrho))
                
            grid.create(par)
            solve.vfi(par)
    
    return

def run_all(par):

    vfi_run(par)
    
    counter = 0
    for i in xrange(0,par.nbeta):
        for j in xrange(0,par.nrho):
                        
            par.het = '' + str(i) + '_' + str(j)
            par.beta = par.beta_mat[i,j]
            par.rho = par.rho_mat[i,j]  

            if par.nbeta*par.nrho > 1:
                counter += 1
                print('\n#het ' + str(counter) + ' of ' + str(par.nbeta*par.nrho))
                
            sim.simulate(par)
    
    return
         
def pref_func(x, par):
    
    if x.size == 4:
        xold = x
        x = np.zeros(5)
        x[0:4] = xold
        x[4] = -1000
    
    par.gbeta = np.zeros(par.nbeta)
    for i in xrange(0,par.nbeta):
        if i == 0:
            par.gbeta[0] = (1.0/par.nbeta)/2.0
        else:
            par.gbeta[i] = par.gbeta[i-1] + (1.0/par.nbeta)

    par.grho = np.zeros(par.nrho)
    for i in xrange(0,par.nrho):
        if i == 0:
            par.grho[0] = (1.0/par.nbeta)/2.0
        else:
            par.grho[i] = par.grho[i-1] + (1.0/par.nrho)
                
    par.beta_mat = np.zeros((par.gbeta.size, par.grho.size))
    par.rho_mat  = np.zeros((par.gbeta.size, par.grho.size))   
    for i, gbeta in enumerate(par.gbeta):
        for j, grho in enumerate(par.grho):
        
            par.beta_mat[i,j] = (0.80 + 0.2*(np.exp(x[0]+np.exp(x[1])*gbeta)) / (1+np.exp(x[0]+np.exp(x[1])*gbeta)))**(0.25)
            par.rho_mat[i,j] = 1 + 10*np.exp(x[2]+np.exp(x[3])*grho+np.exp(x[4])*gbeta) / (1+np.exp(x[2]+np.exp(x[3])*grho+np.exp(x[4])*gbeta))

    return
            
            
def obj_func(x, par):
    
    evals[0] += 1
    if evals[0] > 1:
        print('\nNEW GUESS ' + str(evals[0]))
        print('x')    
        print(x)

    # 1. determine preference matrices
    pref_func(x, par)
    
    print('beta')    
    print(par.beta_mat)
    print('rho')    
    print(par.rho_mat)
    
    # 2. solve and simulate 
    run_all(par)

    # 3. load data  
    t1 = time.time()
            
    sim = plot_sim.SimStruct()        
    for x in ['pop_d', 'pop_n', 'pop_a',  'pop_P', 'pop_Y']:
              
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
    
    # 4. find moments
    t1 = time.time()
    
    t0 = 50
    plot_sim.calc_sim(sim)    
    val = plot_sim.moments(par,sim,t0)
    
    # 5. extra information
    print('puzzle = ' + str(par.sim_puzzle_share) + '')   
    print('borrower = ' + str(par.sim_borrower_share) + '')   
    print('saver = ' + str(par.sim_saver_share) + '')   
    print('corner = ' + str(par.sim_corner_share) + '')  
    
    # return
    print('val = ' + str(val) + '\n')
    return val


def estimate_func(par):
    
    print('\nESTIMATION INITIATED')
        
    print('initial guess:')
    print(par.x)
    
    res = minimize(obj_func, par.x, args=(par), 
                   method='Nelder-Mead', 
                   options={'xtol':1e-2, 'ftol':1e-2, 'disp':False})

    par.x = res.x
    pref_func(res.x, par)
    
    print('ESTIMATION RESULTS')
    print('x')
    print(par.x)
    print('beta')    
    print(par.beta_mat)
    print('rho')    
    print(par.rho_mat)
                            
# main #
if __name__ == '__main__':
    
    # 1. fundamental settings
    case_base = 'Baseline'

    Nhigh = 6 # number of fine robustness checks
    Nlow = 4 # number of rough robustness checks

    d_data = np.array([16.71994019, 0, 0, 0, 0, 13.69159603, 33.11045456, 85.11367798])/100
    a_data = np.array([102.4033051, 0, 0.594442666, 3.388323069, 16.94366837, 59.57957458, 112.6880341, 333.2748413])/100
    n_data = np.array([85.68336487, -60.0524292, -10.41289711, 0, 7.430531502, 48.62488556, 101.6309052, 324.0131226])/100
    
    x_estimate = [-2.5, 1.0, -6.0, 2.0]    
    par = eval('settings.ParStruct_' + case_base + '()') 
    par.nbeta = 5
    par.nrho = 5
    pref_func(np.array(x_estimate), par)
    beta_med = par.beta_mat[2,2]
    rho_med = par.rho_mat[2,2]
    
    # 2. version settings

    versions = [
    
                [2, 'baseline', ['']],

                [1, 'baseline_chi_lose', [1.0, 2.0, 6.0]],
                [1, 'baseline_pi_x_lose_tot', [0.005, 0.01, 0.04]],
                [1, 'baseline_pi_x_gain_tot', [0.03, 0.12, 0.24]],
                [1, 'baseline_u_ast', [0.04, 0.05, 0.06]],
                
                [1, 'baseline_med', ['']],
                [1, 'baseline_med_detail', ['']],

                [1, 'chi_lose', np.linspace(1.00, 8.00, Nhigh)],
                [1, 'pi_x_lose_tot',[0.001, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04]],
                [1, 'pi_x_gain_tot', np.linspace(0.01, 0.95, Nhigh)],

                [1, 'drd', (1+np.linspace(0.4, 0.16, Nlow))**(1/4)-1],          
                [1, 'lamb', np.linspace(0, 0.99, Nhigh)],
                [1, 'varphi', np.linspace(0.25, 1.3, Nlow)],
                [1, 'eta', np.linspace(0, 1, Nlow)],
                                                                        
                [1, 'ra', (1+np.linspace(-0.04, 0.04, Nlow))**(1/4)-1],
                [1, 'Gamma', (1+np.linspace(0.0, 0.04, Nlow))**(1/4)],
              
                [1, 'sigma_psi', np.sqrt(np.linspace(0.0001, 0.08**2, Nlow))],
                [1, 'pi_uu', np.linspace(0.10, 0.50, Nlow)],
                [1, 'mu', np.linspace(0.10, 0.50, Nlow)],             
                [1, 'trans', np.linspace(0.40, 1.1, Nlow)],

                [1, 'sigma_xi', np.sqrt(np.linspace(0.0001, 0.30**2, Nlow))],
                [1, 'u_ast', np.linspace(0.01, 0.14, Nlow)]
                                 
                ]         
                
    # 3. do all                                                                        
    i = 0    
    for version in versions:
                
        evals = np.zeros((1,1))
        do = version[0]
        name = version[1]
        values = version[2]
                        
        for j, value in enumerate(values):

            t1 = time.time()
                
            # a. setup    
            par = eval('settings.ParStruct_' + case_base + '()') 

            par.d_data = d_data            
            par.a_data = a_data
            par.n_data = n_data
            
            par.name = name + '_' + str(j)
            print('\nFOLDER: ' + os.getcwd() + '\\' + par.name)
                                
            # b. specific parameter changes
            par.nbeta = 1
            par.nrho = 1
            par.beta_mat[0,0] = beta_med
            par.rho_mat[0,0] = rho_med
            
            if name in ['baseline']:
                
                par.tableline = 'baseline'
                
                par.nbeta  = 5
                par.nrho   = 5        
                par.x = x_estimate
                pref_func(np.array(par.x), par)
                
                if do != 2:
                    print(par.beta_mat)
                    print(par.rho_mat)
                
            elif name in ['baseline_chi_lose',
                          'baseline_pi_x_lose_tot',
                          'baseline_pi_x_gain_tot',
                          'baseline_u_ast']:
                
                namestriped = name[9:]
                par.tableline = 'par.' + namestriped + ' = ' + str(value)        
                exec(par.tableline)

                par.nbeta  = 5
                par.nrho   = 5        
                par.x = x_estimate
                pref_func(np.array(par.x), par) 
                                
            elif name in ['baseline_med']:
                
                par.tableline = 'baseline_med'
                
            elif name in ['baseline_med_detail']:

                par.T = 160
                par.Ndb = 140
                par.Nnb = 140
                par.Nn = 50       
                par.tableline = 'baseline_med_detail'                
        
            elif name in ['trans']:
           
                par.u_ast = value*par.u_ast
                par.sigma_xi = value*par.sigma_xi
                par.trans = value
                par.tableline = 'par.' + name + ' = ' + str(value)        
                 
            elif name in ['example']:
           
                par.beta_mat[0,0] = 0.92**(1/4)
                par.rho_mat[0,0] = 3.0            
                par.tableline = 'example'
                              
            else:
                
                par.tableline = 'par.' + name + ' = ' + str(value)        
                exec(par.tableline)
                                                
            print("CASE: " + par.tableline)
            
            # d. solve or load
            if do == 1:                
                run_all(par)    
            elif do == 2:
                estimate_func(par)  
                pickle.dump(par, open( os.getcwd() + '\\' + par.name  + '\\Data\\par.p', 'wb' ))
            else:            
                par = pickle.load(open( os.getcwd() + '\\' + par.name + '\\Data\\par.p', 'rb' ))
    
            if name in ['baseline']:
                print(par.x)
                print(par.beta_mat)
                print(par.rho_mat)
                beta_med = par.beta_mat[2,2]
                rho_med = par.rho_mat[2,2]
                x_estimate = par.x
        
            if name in ['baseline_med']:
                grid.allfigs(par)
                plot_sol.allfigs(par)
                
            # e. simulation figures
            if name in ['baseline']:
                plot_sim.allfigs(par,False)
            else:
                plot_sim.allfigs(par)
                
                
            i = i + 1
            del par

    plot_welfare.allfigs()         
    plot_robustness.allfigs()
    plot_sim.print_robustness_table()