import numpy as np
import math
import matplotlib.pyplot as plt
from consav.grids import nonlinspace
import seaborn as sns
sns.set_style("whitegrid")
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]

from consav import linear_interp

# local modules

#############
# lifecycle #
#############

def lifecycle(model,quantiles:bool=False):
    '''
    Plot the lifecycle of the model.
    Keyword arguments:
    quantiles -- if True, plot quantiles instead of mean + quantiles
    '''
    # a. unpack
    par = model.par
    sim = model.sim

    # b. figure
    fig = plt.figure(figsize=(12,12))

    simvarlist = [('p','$p_t$'),
                  ('y','$y_t$'),
                  ('n','$n_t$'),
                  ('m','$m_t$'),
                  ('c','$c_t$'),
                  ('a','$a_t$'),
                  ('discrete','adjuster share')                  
                  ]

    # determine number of rows in figure, given the number of columns
    cols = 2
    rows = math.ceil(len(simvarlist) / cols)

    # x-axis labels
    age = np.arange(par.T)+par.Tmin

    for i,(simvar,simvarlatex) in enumerate(simvarlist):

        ax = fig.add_subplot(rows,cols,i+1)

        simdata = getattr(sim,simvar)[:par.T,:]
        
        # plot
        if quantiles:
            if simvar not in ['discrete','mpc']:
                series = np.percentile(simdata, np.arange(0, 100, 25),axis=1)
                ax.plot(age, series.T,lw=2)
                if i == 0: ax.legend(np.arange(0, 100, 25),title='Quantiles',fontsize=8)
            else:
                ax.plot(age,np.mean(simdata,axis=1),lw=2)

        else:
            ax.plot(age,np.mean(simdata,axis=1),lw=2)
            if simvar not in ['discrete','mpc']:
                ax.plot(age,np.percentile(simdata,25,axis=1),
                    ls='--',lw=1,color='black')
                ax.plot(age,np.percentile(simdata,75,axis=1),
                    ls='--',lw=1,color='black')
        ax.set_title(simvarlatex)
        if par.T > 10:
            ax.xaxis.set_ticks(age[::5])
        else:
            ax.xaxis.set_ticks(age)

        ax.grid(True)
        if i in [len(simvarlist)-i-1 for i in range(cols)]:
            ax.set_xlabel('age')
    plt.savefig('output/life_cycle.png')
    plt.show()

def lifecycle_compare(model1,latex1,model2,latex2,do_euler_errors=False):

    # a. unpack
    par = model1.par
    sim1 = model1.sim
    sim2 = model2.sim

    # b. figure
    fig = plt.figure(figsize=(12,16))

    simvarlist = [('p','$p_t$',None),
                ('n','$n_t$',None),
                ('m','$m_t$',None),
                ('c','$c_t$',None),
                ('d','$d_t$',None),
                ('a','$a_t$',None),
                ('discrete','adjuster share',None)]
            
    if do_euler_errors:
        simvarlist.append(('euler_error_rel','avg. euler error',None))

    age = np.arange(par.T)+par.Tmin
    for i,(simvar,simvarlatex,j) in enumerate(simvarlist):

        ax = fig.add_subplot(4,2,i+1)

        if simvar == 'euler_error_rel':

            simdata = getattr(sim1,simvar)[:par.T-1,:]
            ax.plot(age[:-1],np.nanmean(simdata,axis=1),lw=2,label=latex1)

            simdata = getattr(sim2,simvar)[:par.T-1,:]
            ax.plot(age[:-1],np.nanmean(simdata,axis=1),lw=2,label=latex2)

        else:
            
            simdata = getattr(sim1,simvar)[:par.T,:]
            ax.plot(age,np.mean(simdata,axis=1),lw=2,label=latex1)
            
            simdata = getattr(sim2,simvar)[:par.T,:]
            ax.plot(age,np.mean(simdata,axis=1),lw=2,label=latex2)

        ax.set_title(simvarlatex)
        if par.T > 10:
            ax.xaxis.set_ticks(age[::5])
        else:
            ax.xaxis.set_ticks(age)

        ax.grid(True)
        if simvar in ['discrete','euler_error_rel']:
            if simvar == 'discrete' and not j == 3:
                continue
            ax.set_xlabel('age')
    
        ax.legend()
    
    plt.tight_layout()
    plt.show()

def mpc_over_cash_on_hand(model):
    '''plot mpc as a function of cash-on-hand for given t'''
    p_bar = np.mean(model.sim.p,axis=1)
    n_bar = np.mean(model.sim.n,axis=1)

    c0 = np.zeros(shape=(model.par.T, len(model.par.grid_m)))
    c1 = np.zeros(shape=(model.par.T, len(model.par.grid_m)))
    mpc = np.zeros(shape=(model.par.T, len(model.par.grid_m)))

    m_grid =  nonlinspace(0,model.par.m_max,model.par.Nm,1.1) 

    for t in range(model.par.T):
        t = int(t)    
        for i,m in enumerate(m_grid):
            c0[t,i] = linear_interp.interp_3d(
                    model.par.grid_p,model.par.grid_n,model.par.grid_m,model.sol.c_keep[t],
                    p_bar[t],n_bar[t],m)
            c1[t,i] = linear_interp.interp_3d(
                    model.par.grid_p,model.par.grid_n,model.par.grid_m,model.sol.c_keep[t],  
                    p_bar[t],n_bar[t],m+model.par.mpc_eps)
            mpc[t,i] = (c1[t,i]-c0[t,i])/model.par.mpc_eps

    plt.figure(figsize=(12,8))
    for t in np.arange(5,model.par.T,10):
       plt.plot(model.par.grid_m,np.mean(mpc[t:t+9,:],axis=0),label='t={}-{}'.format(t+model.par.Tmin,t+model.par.Tmin+9))

    plt.xlim(0,5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Cash-on-hand, $m_t$',fontsize=15)
    plt.ylabel('$\mathcal{MPC}_t$',fontsize=15)
    plt.legend(fontsize=15)
    plt.savefig('output/mpc_over_wealth_twoasset.png')
    plt.show()