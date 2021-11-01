import numpy
import matplotlib.pyplot as plt
import random
#import ipdb
import scipy
import time
import copy
import os
import re
import sys
from numpy import loadtxt
from matplotlib import rc
from pylab import rcParams

# the axes attributes need to be set before the call to subplot
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}, size=14)
rc('text', usetex=True)
rc('axes', edgecolor='black', linewidth=0.5)
rc('legend', frameon=False)
rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}'] # \boldmath

#-----------------------------load data-------------------------------------------------------------#
p=300
alpha=0.032
beta=0.010
gamma=0.005
tau_pre=5
tau_post=5
U=0.5
Z=500
kappa=0.25
L=10
sett=0
a=0.05
T=10000

directory='output_p%d_L%d_alpha%.4f_beta%.4f_gamma%.4f_taupre%.1f_kappa%.3f'%(p,L,alpha,beta,gamma,tau_pre,kappa)
print(directory)
if os.path.isdir(directory) == False:
	print("No such directory!")
	exit(1)

t = loadtxt("%s/input_set%d"%(directory,sett))[:,0]
stimulus_id = loadtxt("%s/input_set%d"%(directory,sett))[:,1]
stimulus_amp = loadtxt("%s/input_set%d"%(directory,sett))[:,2]
back_amp = loadtxt("%s/input_set%d"%(directory,sett))[:,3]


indices_L = loadtxt("%s/indicesL_set%d"%(directory,sett))
mmean = loadtxt("%s/mall_set%d"%(directory,sett))[:,1:]

fig, axs = plt.subplots(2, 1, figsize=(6,1.5))

cm = plt.get_cmap('gist_rainbow')
mycolors=[cm(1.*i/(p)) for i in range(p)]
random.Random(1987).shuffle(mycolors)

mycolors=['#E6194B','#4363d8','#ffe119','#3cb44b','#f58231','#911eb4','#42d4f4','#f032e6','#bfef45','#469990']

i=0
for mu in range(L):
	nu=int(indices_L[mu])
	axs[0].plot(t[numpy.where(stimulus_id==nu)], stimulus_amp[numpy.where(stimulus_id==nu)], color=mycolors[i], lw=1)
	i=i+1
#print(back_amp)
axs[0].plot(t, back_amp, color='lightgray')
axs[0].set_ylabel('Input') 

i=0
for mu in range(L):
	nu=int(indices_L[mu])
	axs[1].plot(t, mmean[:,nu], color=mycolors[i], lw=1)
	i+=1
i=0
for mu in range(p):
	if (mu not in indices_L):
		axs[1].plot(t, mmean[:,mu], ls='--', lw=1)
		i=i+1

axs[0].set_xticklabels([])
axs[0].set_xticks([])

axs[0].set_xlim(0,50)
axs[1].set_xlim(0,50)
axs[0].set_ylim(0,1.1)
axs[1].set_ylim(0,1.1)
axs[1].set_ylabel('Overlap') 
axs[1].set_xlabel('Time [sec]')	

fig.savefig("m_ft_set%d_p%d_L%d_alpha%.3f_beta%.3f_gamma%.3f_taupre%.1f_kappa%.3f.png"%(sett,p,L,alpha,beta,gamma,tau_pre,kappa), bbox_inches='tight')