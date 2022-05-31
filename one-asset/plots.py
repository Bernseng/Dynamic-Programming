import numpy as np
from consav.grids import nonlinspace
import tools
import math
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

def plot_consumption(sol,par):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(1,1,1)
    for age in [26, 35, 45, 55, 65, 75, par.T+par.Tmin-1,par.T+par.Tmin] :
        ax.plot(sol.m[age-par.Tmin-1,:],sol.c[age-par.Tmin-1,:], label=f'age = {age}')
    ax.set_xlabel(f"$m_t$")
    ax.set_ylabel(f"$c(m_t)$")
    ax.set_xlim([np.min(par.a_min), 5])
    ax.set_ylim([0,5])
    ax.set_title(f'Consumption function')
    plt.legend()
    plt.show()

def plot_avg_income(sim,par):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(par.T)+par.Tmin+1,np.mean(sim.Y,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Income $Y_t$")
    ax.set_title(f'Average income')
    plt.show()

def plot_avg_cash_on_hand(sim,par):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(par.T)+par.Tmin+1,np.mean(sim.M,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Cash-on-hand $M_t$")
    ax.set_title(f'Average Cash on hands')
    plt.show()

def plot_avg_consumption(sim,par):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(par.T)+par.Tmin+1,np.mean(sim.C,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Consumption $C_t$")
    ax.set_title(f'Average consumption')
    plt.show()

def plot_avg_assets(sim,par):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(par.T)+par.Tmin+1,np.mean(sim.A,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Asset $A_t$")
    ax.set_title(f'Average Asset')
    plt.show()

def lifecycle(par, sim,deciles:bool=False):
    '''
    Plot the lifecycle of the model.
    Keyword arguments:
    deciles -- if True, plot deciles instead of mean + quantiles
    '''

    fig = plt.figure(figsize=(12,12))

    simvarlist = [('P','$P_t$'),
                  ('Y','$Y_t$'),
                  ('M','$M_t$'),
                  ('C','$C_t$'),
                  ('A','$A_t$'),
                  ('mpc','$\mathcal{MPC}_t$'),                  
                  ]

    # determine number of rows in figure, given the number of columns
    cols = 2
    rows = math.ceil(len(simvarlist) / cols)

    # x-axis labels
    age = np.arange(par.T)

    for i,(simvar,simvarlatex) in enumerate(simvarlist):

        ax = fig.add_subplot(rows,cols,i+1)

        simdata = getattr(sim,simvar)[:par.T,:]

        # plot
        if deciles:
            if simvar not in ['discrete','mpc']:
                series = np.percentile(simdata, np.arange(0, 100, 10),axis=1)
                ax.plot(age, series.T,lw=2)
                if i == 0: ax.legend(np.arange(0, 100, 10),title='Deciles',fontsize=8)
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

    plt.show()

def mpc_over_cash_on_hand(model):
    # plot mpc as a function of cash-on-hand for given t

    c0 = np.zeros(shape=(model.par.T, len(model.par.grid_m)))
    c1 = np.zeros(shape=(model.par.T, len(model.par.grid_m)))
    mpc = np.zeros(shape=(model.par.T, len(model.par.grid_m)))

    m_grid =  nonlinspace(0,model.par.m_max,model.par.Nm,1.1) # model.par.grid_m

    for t in range(model.par.T):
        t = int(t)    
        c0[t,:] = tools.interp_linear_1d(m_grid,model.sol.c[t],m_grid[t])
        c1[t,:] = tools.interp_linear_1d(m_grid,model.sol.c[t],m_grid[t]+model.par.mpc_eps)
        for m in m_grid:
            m = int(m)
            mpc[t,m] = (c1[t,m]-c0[t,m])/model.par.mpc_eps

    plt.figure(figsize=(12,8))
    for t in np.arange(0,model.par.T,5):
        plt.plot(model.sol.m[t-model.par.Tmin-1,:],mpc[t,:],linestyle='-',marker='o',label='t={}'.format(t+model.par.Tmin))
    plt.xlim(0,1)
    plt.xlabel('Cash-on-hand, $m_t$')
    plt.ylabel('$\mathcal{MPC}_t$')
    plt.title('$\mathcal{MPC}$ as a function of cash-on-hand (Keep-problem), for mean $p_t$ and $n_t$', fontweight='bold')
    plt.legend()
    plt.show()