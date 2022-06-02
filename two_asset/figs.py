from cProfile import label
from click import style
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from consav.grids import nonlinspace
import seaborn as sns
sns.set_style("whitegrid")
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]
import ipywidgets as widgets

from consav import linear_interp

# local modules
import utility

######################
# decision functions #
######################

def _decision_functions(model,t,i_p,name):

    if name == 'discrete':
        _discrete(model,t,i_p)
    elif name == 'adj':
        _adj(model,t,i_p)
    elif name == 'keep':
        _keep(model,t,i_p)
    elif name == 'post_decision' and t <= model.par.T-2:
        _w(model,t,i_p)        

def decision_functions(model):
    widgets.interact(_decision_functions,
        model=widgets.fixed(model),
        t=widgets.Dropdown(description='t', 
            options=list(range(model.par.T)), value=0),
        i_p=widgets.Dropdown(description='ip', 
            options=list(range(model.par.Np)), value=np.int(model.par.Np/2)),
        name=widgets.Dropdown(description='name', 
            options=['discrete','adj','keep','post_decision'], value='discrete')
        )

def _discrete(model,t,i_p):

    par = model.par

    # a. interpolation
    n, m = np.meshgrid(par.grid_n,par.grid_m,indexing='ij')
    x = m + (1-par.tau)*n
    
    inv_v_adj = np.zeros(x.size)
    linear_interp.interp_1d_vec(par.grid_x,model.sol.inv_v_adj[t,i_p,:,],x.ravel(),inv_v_adj)
    inv_v_adj = inv_v_adj.reshape(x.shape)

    # f. best discrete choice
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1)

    I = inv_v_adj > model.sol.inv_v_keep[t,i_p,:,:]

    x = m[I].ravel()
    y = n[I].ravel()
    ax.scatter(x,y,s=2,label='adjust')
    
    x = m[~I].ravel()
    y = n[~I].ravel()
    ax.scatter(x,y,s=2,label='keep')
        
    ax.set_title(f'optimal discrete choice ($t = {t}$, $p = {par.grid_p[i_p]:.2f}$)',pad=10)

    legend = ax.legend(loc='upper center', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    # g. details
    ax.grid(True)
    ax.set_xlabel('$m_t$')
    ax.set_xlim([par.grid_m[0],par.grid_m[-1]])
    ax.set_ylabel('$n_t$')
    ax.set_ylim([par.grid_n[0],par.grid_n[-1]])
    
    plt.show()

def _adj(model,t,i_p):

    # a. unpack
    par = model.par
    sol = model.sol

    # b. figure
    fig = plt.figure(figsize=(12,6))
    ax_b = fig.add_subplot(1,2,1)
    ax_v = fig.add_subplot(1,2,2)
    
    # c. plot consumption
    ax_b.plot(par.grid_x,sol.d_adj[t,i_p,:],lw=2)
    ax_b.set_title(f'$d^{{adj}}$ ($t = {t}$, $p = {par.grid_p[i_p]:.2f}$)',pad=10)

    # d. plot value function
    ax_v.plot(par.grid_x,sol.inv_v_adj[t,i_p,:],lw=2)
    ax_v.set_title(f'neg. inverse $v^{{adj}}$ ($t = {t}$, $p = {par.grid_p[i_p]:.2f}$)',pad=10)

    # e. details
    for ax in [ax_b,ax_v]:
        ax.grid(True)
        ax.set_xlabel('$x_t$')
        ax.set_xlim([par.grid_x[0],par.grid_x[-1]])

    plt.show()

def _w(model,t,i_p):

    # a. unpack
    par = model.par
    sol = model.sol

    # b. figure
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1,projection='3d')

    n,a = np.meshgrid(par.grid_n, par.grid_a,indexing='ij')

    # c. plot consumption
    ax.plot_surface(n,a,sol.inv_w[t,i_p,:,:],cmap=cm.viridis,edgecolor='none')
    ax.set_title(f'neg. inverse  $w$ ($t = {t}$, $p = {par.grid_p[i_p]:.2f}$)',pad=10)

    # d. details
    ax.grid(True)
    ax.set_xlabel('$n_t$')
    ax.set_xlim([par.grid_n[0],par.grid_n[-1]])
    ax.set_ylabel('$a_t$')
    ax.set_ylim([par.grid_a[0],par.grid_a[-1]])
    ax.invert_xaxis()

    plt.show()

def _keep(model,t,i_p):

    # unpack
    par = model.par
    sol = model.sol

    # b. figure
    fig = plt.figure(figsize=(12,6))
    ax_c = fig.add_subplot(1,2,1,projection='3d')
    ax_v = fig.add_subplot(1,2,2,projection='3d')

    n,m = np.meshgrid(par.grid_n, par.grid_m,indexing='ij')

    # c. plot consumption
    ax_c.plot_surface(n,m,sol.c_keep[t,i_p,:,:],cmap=cm.viridis,edgecolor='none')
    ax_c.set_title(f'$c^{{keep}}$ ($t = {t}$, $p = {par.grid_p[i_p]:.2f}$)',pad=10)

    # d. plot value function
    ax_v.plot_surface(n,m,sol.inv_v_keep[t,i_p,:,:],cmap=cm.viridis,edgecolor='none')
    ax_v.set_title(f'neg. inverse $v^{{keep}}$ ($t = {t}$, $p = {par.grid_p[i_p]:.2f}$)',pad=10)

    # e. details
    for ax in [ax_c,ax_v]:

        ax.grid(True)
        ax.set_xlabel('$n_t$')
        ax.set_xlim([par.grid_n[0],par.grid_n[-1]])
        ax.set_ylabel('$m_t$')
        ax.set_ylim([par.grid_m[0],par.grid_m[-1]])
        ax.invert_xaxis()

    plt.show()

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
    plt.tight_layout()
    plt.savefig('output/life_cycle_twoasset.png')
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
    plt.xlabel('Cash-on-hand, $m_t$')
    plt.ylabel('$\mathcal{MPC}_t$')
    plt.title('$\mathcal{MPC}$ as a function of cash-on-hand (Keep-problem)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('output/mpc_over_wealth_twoasset.png')
    plt.show()