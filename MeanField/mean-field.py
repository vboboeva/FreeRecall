#!/usr/bin/env python
# coding: utf-8


import matplotlib.cm as cm
import os
import tempfile
import numpy as np
import scipy as sp
from scipy.fft import fft, fftfreq
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import time
from tqdm import tqdm

import matplotlib.pyplot as plt
from matplotlib import rc
from pylab import rcParams

from matplotlib import use
use('Agg')

# # the axes attributes need to be set before the call to subplot
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']}, size=14)
# rc('text', usetex=True)
# rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}'] # \boldmath
rc('axes', edgecolor='black', linewidth=0.5)
rc('legend', frameon=False)
rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'


# ======================================================================
# ==================== Dynamical system definitions ====================

def sigma (y):
    out = np.zeros(y.shape)
    out[np.where(y>0)] = 1
    mask = np.where(np.abs(y)<600)
    out[mask] = 1./(1. + np.exp(-y[mask]))
    return out

def rhs (x, beta, gamma, kappa, noise, eps=.1, tau=2.):
    m, n = x
    delta = (2*np.random.rand() - 1.)*noise # np.sqrt(2.*noise)*np.random.normal(size=m.shape) 
    fm = 1./eps * (
        sigma( beta*( m - n - kappa*np.sum(m)**2 + gamma + delta) )
        - m
    )
    fn = 1./tau * ( m - n )

    return fm, fn


# ======================================================================
# ==================== Numerical simulations and plots =================

def simulation(beta, gamma, kappa, noise, T=100, dt=.01, x0=(np.ones(1),np.zeros(1)), eps=.1, tau=2., a=.05, seed=False):
    if seed:
        np.random.seed(seed)
    
    L = len(x0[0])
    
    ts = []
    sol = []
    m,n = tuple(x0)
    for t in range(int(T/dt)):
        ts.append(t*dt)
        sol.append(np.concatenate((m,n)))
        fm, fn = rhs((m,n), beta, gamma, kappa, noise, eps=eps, tau=tau)
        m += dt*fm
        n += dt*fn
            
    return np.array(ts), np.array(sol).T

def classifier (beta, gamma, kappa, noise, T=500, x0=(np.ones(1),np.zeros(1)), dt=0.01, tau=2., eps=.1, skip=None, debug=False, seed=False):
    
    L = len(x0[0])
    
    # solve ODEs
    ts, sol = simulation(beta, gamma, kappa, noise, T=T, dt=dt, x0=x0, tau=tau, eps=eps, seed=seed)
    if skip:
        ts = ts[int(skip/dt):]
        sol = sol[:,int(skip/dt):]
    
    # mean activity over all patterns
    return np.mean(sol[:L])

def classifier_IPR (beta, gamma, kappa, noise, T=500, x0=(np.ones(1),np.zeros(1)), dt=0.01, tau=2., eps=.1, skip=None, debug=False, seed=False):
    
    L = len(x0[0])
    
    # solve ODEs
    ts, sol = simulation(beta, gamma, kappa, noise, T=T, dt=dt, x0=x0, tau=tau, eps=eps, seed=seed)
    if skip:
        ts = ts[int(skip/dt):]
        sol = sol[:,int(skip/dt):]
    
    prob = sol[:L]/np.sum(sol[:L], axis=0)
    IPR = np.sum(prob**2, axis=0)
    
    # mean inverse participation ratio
    return np.mean(IPR)
        
def plot_traj(beta, gamma, kappa, noise, x0=(np.array([1.,0.]),np.array([0.,0.])), dt=.01, T=100, tau=2., eps=.1,
              seed=False, skip_init=False, xticks=None, yticks=None,
              ax=None, filename=None, legend=False):
    
    L=len(x0[0])
    
    ts, sol = simulation(beta, gamma, kappa, noise, x0=x0, T=T, dt=dt, tau=tau, eps=eps, seed=seed)
    
    prob = sol[:L]/np.sum(sol[:L], axis=0)
    IPR = np.sum(prob**2, axis=0)
    
    if not ax:
        fig, ax = plt.subplots(figsize=(5,2.5))
        
    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
    if yticks:
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)
    
    for m,n,i in zip(sol[:L],sol[L:],range(L)):
        ax.plot(ts, m, ls='-',color='C%d'%i)
#         ax.plot(ts, n, ls='--',color='C%d'%i)
#     ax.plot(ts, IPR, ls=':',color='k',lw=2)
    if legend:
        ax.legend(loc="upper right")
    
    if filename:
        plt.savefig(filename,bbox_inches='tight')


def sample_trajectory(L, beta=100, gamma=0., kappa=1.2, noise=0.,tau=2.,T=500):

    m0 = .5*np.random.rand(L)
    n0 = .5*np.random.rand(L)
    m0[0] = 1.
    
    fig, ax = plt.subplots(figsize=(20,5))
    plot_traj(beta, gamma, kappa, noise, x0=(m0,n0),dt=.01, eps=.1, tau=tau,T=T, ax=ax)

    plt.savefig(f"trajectory_L={L}_beta={beta}_gamma={gamma}_kappa={kappa}_noise={noise}_tau={tau}.svg", bbox_inches="tight")
    plt.close(fig)



# ======================================================================
# ======================== Analytical expressions ======================

def gamma_plus (beta, noise=0.):
    return np.log(beta)/beta + noise

def gamma_minus (beta, noise=0.):
    return -np.log(beta)/beta - noise

def n_min(beta, gamma, kappa=0., noise=0.):
    return gamma + (1 + np.log(beta))/beta

def n_max(beta, gamma, kappa=0., noise=0.):
    return 1 + gamma - (1 + np.log(beta))/beta

def t_ON(beta, gamma, kappa=0., noise=0., tau=2.):
    pars = beta, gamma
    return tau*np.log( ( 1 - 1/beta - n_min(*pars) ) / ( 1 - 1/beta - n_max(*pars) ) ) 

def t_OFF(beta, gamma, kappa=0., noise=0., tau=2.):
    pars = beta, gamma
    return tau*np.log( ( - 1/beta + n_max(*pars,kappa=kappa,noise=noise) ) / ( - 1/beta + n_min(*pars,kappa=kappa,noise=noise) ) ) 

def classifier_analytical(beta, gamma, kappa=0., noise=0.,verbose=False):
    
    pars = beta, gamma
    act = t_ON(*pars,kappa=kappa,noise=noise)
    inact = t_OFF(*pars,kappa=kappa,noise=noise)
    
    out = act/(act + inact)
    mask = np.where(np.isnan(act))
    out[mask] = 1.
    mask = np.where(np.isnan(inact))
    out[mask] = 0.
    
    return out

# ======================================================================
# ============================ Phase diagrams ==========================

def phase_diagram_analytical ():
    '''
    Phase diagram in beta-gamma plane from analytical calculations
    (L=1, kappa=0, eps=0)
    '''
    gammas = np.linspace(-.2, .2, 300)
    betas = 10.**np.linspace(1, 3, 300)
    bb, gg = np.meshgrid(betas, gammas)
    phase_an = classifier_analytical(bb, gg)
            
    fig, ax = plt.subplots(figsize=(3,3))
    ax.set_xscale("log")
    _xticks = [10,100,1000]
    _yticks = [-0.2, 0.0, 0.2]
    ax.set_xticks(_xticks)
    ax.set_yticks(_yticks)
    ax.set_yticklabels(_yticks)
    ax.tick_params(top=True,bottom=True,left=True,right=True,which='both',pad=10)
    im = ax.contourf(bb, gg, phase_an, np.linspace(0,1,300), cmap='jet')#, vmin=0, vmax=1)

    upper = gamma_plus(betas)
    lower = gamma_minus(betas)
    mask = np.where(upper < 0.2)
    ax.plot(betas[mask],upper[mask],ls='--',color='white', linewidth=3)
    ax.plot(betas[mask],lower[mask],ls='--',color='white', linewidth=3)
    plt.savefig('phase_analytical.svg',bbox_inches="tight")


def phase_diagrams_numerical ():

    gammas = np.linspace(-.2, .2, 50)
    betas = 10.**np.linspace(1, 3, 50)

    Ls = [1, 1, 1, 1, 2, 2]
    kappas = [0, .05, 0, .05, .05, .05]
    noises = [0, 0, .1, .1, 0, .1]

    for kappa,noise,L in zip(kappas,noises,Ls):
        print("L =", L, "  noise =", noise, "  kappa =", kappa)
        filename = 'mean_activity_L=%d_kappa=%.2f_noise=%.2f'%(L,kappa,noise)

        def func(p):
            m0 = 1 - .05*np.random.rand(L)
            n0 = .05*np.random.rand(L)
            m0[0] = 1.
            return classifier(*p,T=500, x0=(m0,n0))

        try:
            '''
            Load phase diagram if already been computed

            '''
            phases = np.load(f"{filename}.npy")

            if phases.shape != (len(gammas),len(betas)):
                raise FileNotFoundError

        except FileNotFoundError:
            '''
            Compute again phase diagram if not saved or if it has the wrong shape

            '''

            phases = np.zeros((len(gammas), len(betas)))
            for g, gamma in enumerate(gammas):
                print(f'  {g}/{len(gammas)}', end='\r')

                args = [(beta, gamma, kappa, noise) for beta in betas]

                results = map(func, args)

                # # with multiprocessing
                # with ProcessPoolExecutor(max_workers=50) as executor:
                #     results = executor.map(func, args)

                phases[g] = np.array([res for res in results])

            np.save(f"{filename}.npy", phases)

        bb, gg = np.meshgrid(betas, gammas)
        fig, ax = plt.subplots(figsize=(3,3))
        ax.set_xscale('log')
        im = ax.contourf(bb, gg, phases, np.linspace(0,1,300), vmin=0, vmax=1, cmap='jet')
        plt.savefig(f'{filename}.svg')



# ======================================================================
# ====================== Dynamics at selected points ===================

def sample_dynamics_selected_point ():

    gmax = 0.105
    gmin = -gmax
    beta = 400
    kmax = 0.05
    nmax = 0.1

    points = [
        # Figure S1-B
        (1, (beta, gmin, 0, 0), 100),    # paramagnetics
        (1, (beta, 0,    0, 0), 100),    # deterministic oscillations
        (1, (beta, gmax, 0, 0), 100),    # ferromagnetic
        # 
        (1, (beta, gmin, 0, nmax), 100), # random oscillations
        (1, (beta, 0,    0, nmax), 100), # random activations
        (1, (beta, gmax, 0, nmax), 100), # random inactivations
        #
        (2, (beta, gmin, kmax, 0), 100), # paramagnetic
        (2, (beta, 0,    kmax, 0), 100), # latching
        (2, (beta, gmax, kmax, 0), 100), # latching
        #
        # Figure S3
        (2, (100, 0, 0, 0), 150),        # independent
        (2, (100, 0, kmax, 0), 150),     # latching - out-of-phase
        (16, (100, 0, 10*kmax, 0), 50),  # latching
        (16, (100, 0, 10*kmax, 1), 50),  # distraction
    ]

    _xticks = [0, 25, 50, 75, 100, 125, 150]
    _yticks = [0, 0.5, 1.0]
    
    for p in points:
        fig, ax = plt.subplots()
        L, pars, T = p
        beta, gamma, kappa, noise = pars
        ax.set_xticks(_xticks)
        ax.set_yticks(_yticks)
        ax.set_ylim([-.1,1.1])
        x0 = ( np.ones(L) - 0.001*np.random.rand(L), np.zeros(L) + 0.001*np.random.rand(L) )
        plot_traj(*pars,T=T,legend=False, x0=x0, ax=ax, tau=2., eps=.05)
        plt.savefig(f'sample_L={L}_noise={noise:.3f}_kappa={kappa:.3f}_beta={beta}_gamma={gamma:.3f}.svg', bbox_inches='tight')


if __name__ == "__main__":

    # Figure S2-A
    phase_diagram_analytical()

    # Panels in Figure S1-A and Figure S2-B
    phase_diagrams_numerical()

    # Panels in Figure S1-B and Figure S3
    sample_dynamics_selected_point()
