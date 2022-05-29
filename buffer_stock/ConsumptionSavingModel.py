# -*- coding: utf-8 -*-
"""ConsumptionSavingModel

Solves the Deaton-Carroll buffer-stock consumption model with vfi or egm:

"""

##############
# 1. imports #
##############

import time
import numpy as np
from numba import njit, prange
from scipy import optimize

# consav package
from consav import ModelClass, jit # baseline model class and jit
from consav import linear_interp # linear interpolation
from consav.grids import nonlinspace # grids
from consav.quadrature import log_normal_gauss_hermite # income shocks
from consav.misc import elapsed
from quadrature import create_PT_shocks

import figs
import simulate
import trans
import quadrature
import utility
import last_period
import egm

############
# 2. model #
############

class ConsumptionSavingModelClass(ModelClass):
    
    #########
    # setup #
    #########      

    def settings(self):
        """ fundamental settings """

        # namespaces
        self.namespaces = []
        
        # other attributes
        self.other_attrs = []
        
        # savefolder
        self.savefolder = 'saved'
        
        # for safe type inference
        self.not_floats = ['solmethod','T','TR','Nxi','Npsi','Nm','Na','simN','Nshocks','sim_seed']
    
    def setup(self):
        """ set baseline parameters """

        par = self.par
        
        # horizon and life cycle
        par.Tmin = 25
        par.T = 80 - par.Tmin
        par.Tr = 65 - par.Tmin # retirement age (end-of-period), no retirement if TR = T
        par.G = 1.02
        par.L = np.ones(par.T-1)
        par.L[0:par.Tr] = np.linspace(1,1/par.G,par.Tr) 
        par.L[par.Tr-1] = 0.67 # drop in permanent income at retirement age
        par.L[par.Tr-1:] = par.L[par.Tr-1:]/par.G # constant permanent income after retirement
        
        # preferences
        par.rho = 2.0 # CRRA coeficient
        par.beta = 0.96 # discount factor

        # returns and incomes
        par.R = 1.03
        par.sigma_psi = 0.1
        par.sigma_xi = 0.1
        par.Npsi = 5
        par.Nxi = 5
        par.mpc_eps = 0.00855
        
        # grids
        par.Nm = 100
        par.m_max = 10
        par.m_phi = 1.1 #curvature parameter
        par.Na = 100
        par.a_max = par.m_max+1.0
        par.a_phi = 1.1 #curvature parameter
        par.Np = 50
        par.p_min = 1e-4
        par.p_max = 3.0
        
        # simulation
        par.sigma_m0 = 0.2
        par.mu_m0 = 0.0
        par.mu_p0 = 0.0
        par.sigma_p0 = 0.2
        par.simN = 10000 # number of persons in simulation
        #par.simT = 80 # number of periods in simulation
        par.sim_seed = 1998
        par.euler_cutoff = 0.02
        
        # misc
        par.t = 0
        par.tol = 1e-8
        par.do_print = False
        par.do_print_period = False
        par.do_marg_u = False

    def allocate(self):
        """ allocate model, i.e. create grids and allocate solution and simluation arrays """

        self.create_grids()
        self.solve_prep()
        self.simulate_prep()

    def create_grids(self):
        """ construct grids for states and shocks """

        par = self.par

        # a. post-decision states
        par.grid_a = np.ones((par.T,par.Na))
        par.a_min = np.zeros(par.T) # never any borriwng
        for t in range(par.T):
            par.grid_a[t,:] = nonlinspace(par.a_min[t]+1e-6,par.a_max,par.Na,par.a_phi)
        
        # b. states
        par.grid_m = np.ones((par.T,par.Nm))
        for t in range(par.T):
            par.grid_m[t,:] = nonlinspace(par.a_min[t]+1e-6,par.m_max,par.Nm,par.m_phi)
        
        # c. shocks
        shocks = create_PT_shocks(
            sigma_psi=par.sigma_psi,
            Npsi=par.Npsi,
            sigma_xi=par.sigma_xi,
            Nxi=par.Nxi,
            )
        par.psi,par.psi_w,par.xi,par.xi_w,par.Nshocks = shocks
        par.w = par.psi_w*par.xi_w # weights
        assert 1-np.sum(par.w) < 1e-8 # == summing to 1

        # d. set seed
        np.random.seed(par.sim_seed)
        
    """
        # e. conditions
        par.FHW = np.float(par.G/par.R)
        par.AI = np.float((par.R*par.beta)**(1/par.rho))
        par.GI = np.float(par.AI*np.sum(par.w*par.psi**(-1))/par.G)
        par.RI = np.float(par.AI/par.R)
        #par.WRI = np.float(par.pi**(1/par.rho)*par.AI/par.R)
        par.FVA = np.float(par.beta*np.sum(par.w*(par.G*par.psi)**(1-par.rho)))

        # g. check for existance of solution
        self.check(do_print=False)
    
    def check(self,do_print=True):
        #check parameters 

        par = self.par

        if do_print:
            print(f'FHW = {par.FHW:.3f}, AI = {par.AI:.3f}, GI = {par.GI:.3f}, RI = {par.RI:.3f}, WRI = {par.WRI:.3f}, FVA = {par.FVA:.3f}')

        # check for existance of solution
        if par.sigma_xi == 0 and par.sigma_psi == 0: # no risk
            if par.GI >= 1 and par.RI >= 1:
                raise Exception('GI >= 1 and RI >= 1: no solution')
        else:
            if par.FVA >= 1 or par.WRI >= 1:
                raise Exception('FVA >= 1 or WRI >= 1: no solution')
     """                                
    #########
    # solve #
    #########
    
    def precompile_numba(self):
        """ solve the model with very coarse grids and simulate with very few persons"""

        par = self.par

        tic = time.time()

        # a. define
        fastpar = dict()
        fastpar['do_print'] = False
        fastpar['do_print_period'] = False
        fastpar['T'] = 2
        fastpar['Np'] = 3
        fastpar['Nm'] = 3
        fastpar['Na'] = 3
        fastpar['simN'] = 2

        # b. apply
        for key,val in fastpar.items():
            prev = getattr(par,key)
            setattr(par,key,val)
            fastpar[key] = prev

        self.allocate()

        # c. solve
        self.solve()

        # d. simulate
        self.simulate()

        # e. reiterate
        for key,val in fastpar.items():
            setattr(par,key,val)
        
        self.allocate()

        toc = time.time()
        if par.do_print:
            print(f'numba precompiled in {toc-tic:.1f} secs')

    
    def solve_prep(self):
        """ allocate model, i.e. create grids and allocate solution and simluation arrays """
        
        par = self.par
        sol = self.sol

        sol_shape = (par.T,par.Na+1)
        sol.m = np.zeros(sol_shape)
        sol.c = np.zeros(sol_shape)
        sol.inv_v = np.zeros(sol_shape)

        
    def solve(self,do_print=True):
        """ gateway for solving the model """

        par = self.par

        # a. reset solution and simulation arrays
        #self.allocate()

        # b. solve
        tic = time.time()

        # backwards induction
        for t in reversed(range(self.par.T)):
            self.par.t = t
            
            with jit(self) as model:
                par = model.par
                sol = model.sol
                
                m = np.zeros(par.Na)
                c = np.zeros(par.Na)
                inv_v = np.zeros(par.Na)
                
                # last period
                if t == par.T-1:
                    last_period.last_period(sol,par)
                
                # other periods    
                else:
                    egm.egm(par,sol,t,m,c,inv_v) # solve by egm

                    # ii. add zero consumption
                    sol.m[t,0] = par.a_min[t]
                    sol.m[t,1:] = m
                    sol.c[t,0] = 0
                    sol.c[t,1:] = c
                    sol.inv_v[t,0] = 0
                    sol.inv_v[t,1:] = inv_v
            
        toc = time.time()

        if par.do_print:
            print(f'model solved in {toc-tic:.1f} secs')


    ############
    # simulate #
    ############
    
    def simulate_prep(self):
        
        par = self.par
        sim = self.sim
        
        # initial
        sim.m0 = np.zeros(par.simN)
        sim.p0 = np.zeros(par.simN)
        
        sim.utility = np.zeros(par.simN) # necessary?
        
        # states and choices
        sim_shape = (par.T,par.simN)
        sim.m = np.zeros(sim_shape)
        sim.y = np.zeros(sim_shape)
        sim.p = np.zeros(sim_shape)
        
        sim.P = np.zeros(sim_shape)
        sim.Y = np.zeros(sim_shape)
        sim.M = np.zeros(sim_shape)
        
        sim.c = np.zeros(sim_shape)
        sim.a = np.zeros(sim_shape)
        
        sim.C = np.zeros(sim_shape)
        sim.A = np.zeros(sim_shape)

        # shocks
        sim.psi = np.zeros(sim_shape)
        sim.xi = np.zeros(sim_shape)
        
        # mpc
        sim.mpc = np.zeros(sim_shape)
        
        # euler
        euler_shape = (par.T-1,par.simN)
        sim.euler_error = np.zeros(euler_shape)
        sim.euler_error_c = np.zeros(euler_shape)
        sim.euler_error_rel = np.zeros(euler_shape)
    
    def simulate(self,do_print=True,seed=2017):
        """ simulate the model """

        par = self.par
        sol = self.sol
        sim = self.sim
        
        tic = time.time()
        
        # initial states
        sim.m0[:] = np.random.lognormal(mean=par.mu_m0,sigma=par.sigma_m0,size=par.simN)
        sim.p0[:] = np.random.lognormal(mean=par.mu_p0,sigma=par.sigma_p0,size=par.simN)

        # shocks
        I = np.random.choice(par.Nshocks,
            size=(par.T,par.simN),
            p=par.w)
        sim.psi[:,:] = par.psi[I]
        sim.xi[:,:] = par.xi[I]
        
        # call
        with jit(self) as model:

            par = model.par
            sol = model.sol
            sim = model.sim
            
            simulate.life_cycle(par,sol,sim)

        toc = time.time()

        # e. renormalized
        sim.P[:,:] = np.exp(sim.p)
        sim.Y[:,:] = np.exp(sim.y)
        sim.M[:,:] = sim.m*sim.P
        sim.C[:,:] = sim.c*sim.P
        sim.A[:,:] = sim.a*sim.P
        

        if do_print:
            print(f'model simulated in {toc-tic:.1f} secs')
        
        # euler errors
    
    #########
    # plots #
    #########
    """
    plot_value_function_convergence = figs.plot_value_function_convergence
    plot_consumption_function_convergence = figs.plot_consumption_function_convergence
    plot_consumption_function_convergence_age = figs.plot_consumption_function_convergence_age
    plot_consumption_function_pf = figs.plot_consumption_function_pf
    plot_buffer_stock_target = figs.plot_buffer_stock_target
    plot_simulate_cdf_cash_on_hand = figs.plot_simulate_cdf_cash_on_hand
    plot_simulate_consumption_growth = figs.plot_simulate_consumption_growth
    plot_life_cycle_income = figs.plot_life_cycle_income
    plot_life_cycle_cashonhand = figs.plot_life_cycle_cashonhand
    plot_life_cycle_consumption = figs.plot_life_cycle_consumption
    plot_life_cycle_assets = figs.plot_life_cycle_assets
    """