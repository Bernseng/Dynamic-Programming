import numpy as np
from consav.grids import nonlinspace
import tools
import math
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

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
                  ('A','$A_t$')                  
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
    plt.tight_layout()
    plt.show()
    plt.savefig('output/lifecycle_oneasset.png')
    

def mpc_over_cash_on_hand(par,sol,sim):
    """ plot mpc as a function of cash-on-hand for given t """
    
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
    
    plt.figure(figsize=(9,6))
    for t in np.arange(5,par.T,10):
        plt.plot(m_grid[1:],np.mean(mpc[t:t+9,1:],axis=0),label='t={}-{}'.format(t+par.Tmin,t+par.Tmin+9),lw=2.3)
    plt.xlim(0,3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Cash-on-hand, $m_t$',fontsize=15)
    plt.ylabel('$\mathcal{MPC}_t$',fontsize=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.savefig('output/mpc_over_wealth_oneasset.png')
    plt.show()
    
def mpc_over_cash_on_hand_60(par,sol,sim):
    """ plot mpc as a function of cash-on-hand for 60-69 to check retirement behavior """

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

    plt.figure(figsize=(9,6))
    for t in np.arange(par.T-20,par.T-10,1):
        plt.plot(m_grid[1:],mpc[t,1:],label='t={}'.format(t+par.Tmin),lw=2.3)
    plt.xlim(0,3)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Cash-on-hand, $m_t$',fontsize=15)
    plt.ylabel('$\mathcal{MPC}_t$',fontsize=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.show()
    plt.savefig('output/mpc_over_wealth_oneasset_60s.png')
