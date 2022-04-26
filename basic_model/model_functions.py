import numpy as np
from scipy import optimize
from timeit import default_timer as timer
from utils import gauss_hermite_lognormal

__all__ = ['settings','setup','allocate','solve','simulate']

# settings function
def settings(model):
    """ choose fundamental settings """
    
    # a. namespaces
    model.namespaces = ['par','sol','sim']

# setup function
def setup(model):
    """ set baseline parameters (independent elements) """
        
    par = model.par

    # a. parameters
    par.T = 100 # number of periods
    par.rho = 2.0 # CRRA parameter
    par.beta = 0.96 # discount factor
    par.sigma_xi = 0.1 # std. of transitory shock
    par.Nxi  = 8 # number of quadrature points for xi
    par.R = 1.04 # return factor
    par.Nm = 100 # number of nodes in grid for m
    par.m_min = 1e-6 # minimum point in grid for m
    par.m_max = 20.0 # maximum point in grid for m
    par.m_phi = 1.3 # curvature parameter for grid for m
    par.simN = 10_000 # number of persons in simulation


# allocate space for solutions and simulations
def allocate(model):
    """ allocate model, i.e. create grids and allocate solution and simluation arrays (dependent elements) """
    
    par = model.par
    sol = model.sol
    sim = model.sim   
    
    # a. grids
    # non-equally spaced grid
    par.grid_m = np.linspace(par.m_min,par.m_max,par.Nm)
    
    # nodes and weights using Gauss-Hermite
    par.xi, par.xi_w = gauss_hermite_lognormal(sigma=par.sigma_xi,n=par.Nxi) 
    
    # b. solution
    sol.c = np.zeros((par.T,par.Nm))
    sol.v = np.zeros((par.T,par.Nm))
    
    # c. simulation
    sim.m = np.zeros((par.T,par.simN))
    sim.c = np.zeros((par.T,par.simN))
    sim.a = np.zeros((par.T,par.simN))
    sim.xi = np.zeros((par.T,par.simN))

# solve functions
def value_of_choice(c,t,m,model):
    """ value of choice """
    
    par = model.par
    sol = model.sol
    
    if c <= 0: return -np.inf
    
    # a. end-of-period assets
    a = m - c
    
    # b. next-period cash-on-hand
    m_next = par.R*a + par.xi
    
    # c. continuation value
    v_next = np.interp(m_next,par.grid_m,sol.v[t+1,:],) 
            
    # d. value-of-choice
    total = c**(1-par.rho)/(1-par.rho) + par.beta*np.sum(par.xi_w*v_next)
    return total    
    
def solve(model):
    """ solve model with vfi """
    
    par = model.par
    sol = model.sol
    
    start = timer()
    
    # a. last-period
    sol.c[-1,:] = par.grid_m
    sol.v[-1,:]= par.grid_m**(1-par.rho)/(1-par.rho) 
    
    # b. backwards induction
    for t in reversed(range(par.T-1)):
        for i_m,m in enumerate(par.grid_m):

            obj = lambda c: -value_of_choice(c,t,m,model)
            result = optimize.minimize_scalar(obj,method='bounded',bounds=(0,m))

            sol.c[t,i_m] = result.x
            sol.v[t,i_m]= -result.fun 
    print(f'model solved in {timer()-start:.1f} seconds')

# simulate function
def simulate(model,seed=1917):
    """ simulate model """
    
    par = model.par
    sol = model.sol
    sim = model.sim
    
    start = timer()
    
    # a. shocks
    if not seed is None: np.random.seed(seed)
    sim.xi[:,:] = np.exp(np.random.normal(-0.5*par.sigma_xi**2,par.sigma_xi,size=(par.T,par.simN)))
    
    # a. initial cash-on-hand
    sim.m[0,:] = sim.xi[0,:] 
    
    # c. time loop
    for t in range(par.T):
        
        # i. consumption
        sim.c[t,:] = np.interp(sim.m[t,:],par.grid_m,sol.c[t],)

        # ii. end-of-period assets
        sim.a[t,:] = sim.m[t,:]-sim.c[t,:]

        # iii. next-period
        if t < par.T-1: 
            sim.m[t+1] = par.R*sim.a[t,:] + sim.xi[t+1,:]
            
    print(f'model simulated in {timer()-start:.1f} seconds')    