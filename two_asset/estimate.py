# Import package
import numpy as np
from regex import P
import scipy.optimize as optimize
from consav import linear_interp # for linear interpolation
from TwoAssetModel import TwoAssetModelClass as model

def updatepar(par, parnames, parvals):

    for i,parval in enumerate(parvals):
        parname = parnames[i]
        setattr(par,parname,parval) # It gives the attibute parname the new value parval, within the par class
    return par

def calc_moments(par,data):
    #agegrid = np.arange(par.moments_minage,par.moments_maxage+1)-par.Tmin+1 # define the cell which correspond to the age we want the mean for. e.g. age 40-55 --> agegrid: 16-31

    #return np.array([np.mean(data.a,1),np.mean(data.y,1)]) # both a and y
    return np.mean(data.a,1) # only a

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
    moments = np.nan + np.zeros((data.moments.size,par.moments_numsim)) # both a
    #moments = np.nan + np.zeros((data.moments.shape[0],data.moments.shape[1],par.moments_numsim)) # both a and y

    for s in range(par.moments_numsim):

        # Simulate
        model.simulate()

        # Calculate moments
        moments[:,s] = calc_moments(par,model.sim)    # only a
        #moments[:,:,s] = calc_moments(par,model.sim) * scale  # both a and y

    # Mean of moments         
    moments = np.mean(moments,1)    # only a
    #moments = np.mean(moments,axis=2)   # both a and y

    # Objective function
    if hasattr(par, 'weight_mat'):
        weight_mat_inv = np.linalg.inv(par.weight_mat)  
    else:
        weight_mat_inv = np.eye(moments.size)   # The identity matrix and I^-1=I
    
    diff = (moments-data.moments).reshape(moments.size,1)

    return (diff.T@weight_mat_inv @diff).ravel()
