import numpy as np
import math

def gauss_hermite(n):
    """ computing nodes and weight using Gauss-Hermite quadrature """
    
    # a. calculations
    i = np.arange(1,n)
    a = np.sqrt(i/2)
    CM = np.diag(a,1) + np.diag(a,-1)
    L,V = np.linalg.eig(CM)
    I = L.argsort()
    V = V[:,I].T

    # b. nodes and weights
    x = L[I]
    w = np.sqrt(math.pi)*V[:,0]**2

    return x,w
    
def gauss_hermite_lognormal(sigma,n): 
    
    x, w = gauss_hermite(n)
    
    # transform to log normal shocks
    xlog = np.exp(sigma*np.sqrt(2)*x)
    wlog = w/np.sqrt(np.pi)

    return xlog, wlog