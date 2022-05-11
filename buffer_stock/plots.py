import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

def plot_consumption(sol,par):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(1,1,1)
    for age in [26, 35, 45, 55, 65, 75, par.T+par.age_min-1,par.T+par.age_min] :
        ax.plot(sol.m[age-par.age_min-1,:],sol.c[age-par.age_min-1,:], label=f'age = {age}')
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
    ax.plot(np.arange(par.simT)+par.age_min+1,np.mean(sim.Y,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Income $Y_t$")
    ax.set_title(f'Average income')
    plt.show()

def plot_avg_cash_on_hand(sim,par):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(par.simT)+par.age_min+1,np.mean(sim.M,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Cash-on-hand $M_t$")
    ax.set_title(f'Average Cash on hands')
    plt.show()

def plot_avg_consumption(sim,par):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(par.simT)+par.age_min+1,np.mean(sim.C,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Consumption $C_t$")
    ax.set_title(f'Average consumption')
    plt.show()

def plot_avg_assets(sim,par):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)
    ax.plot(np.arange(par.simT)+par.age_min+1,np.mean(sim.A,1))
    ax.set_xlabel(f"age")
    ax.set_ylabel(f"Asset $A_t$")
    ax.set_title(f'Average Asset')
    plt.show()