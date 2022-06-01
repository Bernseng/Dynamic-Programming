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
                  #('mpc','$\mathcal{MPC}_t$'),                  
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
        if deciles:
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

    plt.show()

def mpc_over_cash_on_hand(par,sol,sim):
    # plot mpc as a function of cash-on-hand for given t
    m_grid =  nonlinspace(0,par.m_max,par.Nm+1,1.1) # par.grid_m

    c0 = np.zeros(shape=(par.T, len(m_grid)))
    c1 = np.zeros(shape=(par.T, len(m_grid)))
    mpc = np.zeros(shape=(par.T, len(m_grid)))

    for t in range(par.T):
        t = int(t)    

        c0[t,:] = tools.interp_linear_1d(sol.m[t,:],
                                        sol.c[t,:],#*np.exp(sim.p[t,:]),
                                        m_grid#*np.exp(sim.p[t,:])
                                        )
        bump = par.mpc_eps / np.mean(sim.P[t,:])
        c1[t,:] = tools.interp_linear_1d(sol.m[t,:],
                                        sol.c[t,:],#*np.exp(sim.p[t,:]),
                                        bump+m_grid,#*np.exp(sim.p[t,:])
                                        )

        for i,m in enumerate(m_grid):
            if i == 0: continue
            mpc[t,i] = (c1[t,i]-c0[t,i])/bump
    
    plt.figure(figsize=(12,8))
    for t in np.arange(0,par.T-10,10):
        plt.plot(m_grid[1:],np.mean(mpc[t:t+9,1:],axis=0),label='t={}-{}'.format(t+par.Tmin,t+par.Tmin+9))
    plt.xlim(0,3)
    plt.xlabel('Cash-on-hand, $m_t$')
    plt.ylabel('$\mathcal{MPC}_t$')
    plt.title('$\mathcal{MPC}$ as a function of cash-on-hand', fontweight='bold')
    plt.legend()
    plt.show()