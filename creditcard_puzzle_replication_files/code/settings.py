# -*- coding: utf-8 -*-
"""
Settings module.

version: 1.0.
@author: Jeppe Druedahl, 2017.
"""

# imports #
from __future__ import division, print_function, absolute_import

import os
import numpy as np

# classes #
class ParStruct_Baseline:
    def __init__(self):

        ###########
        # general #
        ###########

        self.name = 'baseline'
                 
        self.T = 120
        self.epsilon = 10**(-8)
                 
                 
        ###########                 
        # utility #
        ###########
        
        self.nbeta = 1;   
        self.nrho = 1;        
        self.beta_mat = np.zeros((self.nbeta,self.nrho))
        self.rho_mat = np.zeros((self.nbeta,self.nrho))        
        
        self.beta_mat[0,0] = 0.90**0.25
        self.rho_mat[0,0] = 3.0
        
        self.pi_death = 0.01


        ##########        
        # income #
        ##########
        
        self.Gamma = 1.02**(1/4)
               
        self.u_ast = 0.07  
        self.mu = 0.30
        self.pi_uu = 0.07
        
        self.sigma_psi = np.sqrt(0.01*4/11)
        self.sigma_xi = np.sqrt(0.01*4)


        ###################
        # assets and debt #
        ###################
        
        self.ra = (1.00  - 0.0148)**(1/4)-1
        self.drd = (1.00 + 0.1236)**(1/4)-1
        
        self.lamb = 0.03
        self.eta = 0.0
        self.varphi = 0.74
        self.chi = 1.0
        
        self.trans = 1.0
        
        self.x_additive = False       

        self.pi_x_lose_tot = 0.0264
        self.pi_x_gain_tot = 0.0607

        self.chi_lose = 4.0
        self.chi_gain = 0.0
                
        self.no_x = 0


        ########        
        # grid #
        ########

        # shocks
        self.Nxi = 4
        self.Npsi = 4

        # states
        self.Nx = 2
        self.Nu = 2

        self.Ndb = 80
        self.Upsilon = 6
        
        self.Nnb = 80
        self.nb_max = 4
        self.phi_nb = 1.1
        self.nb_end = [4.5, 5, 5.5, 6]

        # choices        
        self.Nn = 40            
        self.phi_n = 1.2

        self.d_stepsize = 0.0050


        ############
        # simulate #
        ############
        
        self.simT = 600        
        self.simBurnIn = 500
        self.N = 50000
