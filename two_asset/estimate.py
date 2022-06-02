# Import package
import numpy as np
from regex import P
import scipy.optimize as optimize
from consav import linear_interp # for linear interpolation
from DurableConsumptionModel import DurableConsumptionModelClass as model

def maximum_likelihood(model, est_par, theta0, data,do_stderr):
    
    # Check the parameters
    assert (len(est_par)==len(theta0)), 'Number of parameters and initial values do not match'
    
    #Estimation
    obj_fun = lambda x: -log_likelihood(x,model,est_par,data)
    res = optimize.minimize(obj_fun,theta0)

    return res

def log_likelihood(theta, model, est_par, data):
    
    #Update parameters
    par = model.par
    sol = model.sol

    par = updatepar(par,est_par,theta)

    # Solve the model
    model.create_grids()
    model.solve()

    # Predict consumption
    t = data.t 
    c_predict = linear_interp.linear_interp_1d(sol.m[t,:],sol.c[t,:],data.m)
    C_predict = c_predict*data.P        #Renormalize

    # Calculate errors
    error = data.logC -np.log(C_predict)

    # Calculate log-likelihood
    log_lik_vec = - 0.5*np.log(2*np.pi*par.sigma_eta**2)
    log_lik_vec += (- (error**2)/(2*par.sigma_eta**2) )
    
    return np.mean(log_lik_vec) 

def updatepar(par, parnames, parvals):

    for i,parval in enumerate(parvals):
        parname = parnames[i]
        setattr(par,parname,parval) # It gives the attibute parname the new value parval, within the par class
    return par

def calc_moments(par,data):
    #agegrid = np.arange(par.moments_minage,par.moments_maxage+1)-par.Tmin+1 # define the cell which correspond to the age we want the mean for. e.g. age 40-55 --> agegrid: 16-31
    noise_a = data.a + np.random.normal(0,par.moments_noise,size=data.a.shape)  # introduce noise to the data on top of new realizations of shocks
    noise_y = data.y + np.random.normal(0,par.moments_noise,size=data.y.shape)  # introduce noise to the data on top of new realizations of shocks
    
    #return np.array([np.mean(noise_a[agegrid,:],1),np.mean(noise_y[agegrid,:],1)])
    return np.array([np.mean(noise_a,1),np.mean(noise_y,1)]) # both a and y
    #return np.mean(noise_a[agegrid,:],1) # only a

def method_simulated_moments(model,est_par,theta0,data):

    # Check the parameters
    assert (len(est_par)==len(theta0)), 'Number of parameters and initial values do not match'
    
    # Calculate data moments
    data.moments = calc_moments(model.par,data)

    # Estimate
    obj_fun = lambda x: sum_squared_diff_moments(x,model,est_par,data)
    res = optimize.minimize(obj_fun,theta0, method='BFGS')

    return res

def sum_squared_diff_moments(theta0,model,est_par,data,scale=1):

    par = model.par
    #Update parameters
    par = updatepar(par,est_par,theta0)

    # Solve the model
    model.create_grids()
    model.solve()

    # Simulate the momemnts
    #moments = np.nan + np.zeros((data.moments.size,par.moments_numsim)) # both a
    moments = np.nan + np.zeros((data.moments.shape[0],data.moments.shape[1],par.moments_numsim)) # both a and y

    for s in range(par.moments_numsim):

        # Simulate
        model.simulate()

        # Calculate moments
        #moments[:,s] = calc_moments(par,model.sim)    # only a
        moments[:,:,s] = calc_moments(par,model.sim) * scale  # both a and y

    # Mean of moments         
    #moments = np.mean(moments,1)    # only a
    moments = np.mean(moments,axis=2)   # both a and y

    # Objective function
    if hasattr(par, 'weight_mat'):
        weight_mat_inv = np.linalg.inv(par.weight_mat)  
    else:
        weight_mat_inv = np.eye(moments.size)   # The identity matrix and I^-1=I
    
    diff = (moments-data.moments).reshape(moments.size,1)

    return (diff.T@weight_mat_inv @diff).ravel()
