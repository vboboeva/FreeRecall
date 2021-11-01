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
dJ=0.01
J_EI=1.2

N=5000
a=0.02

p=600
U=0.5
Z=500
L=16
sett=0
T=500000

directory='output_N%d_a%.2f_p%d_L%d_dJ%.4f_JEI%.1f'%(N,a,p,L,dJ,J_EI)

print(directory)

if os.path.isdir(directory) == False:
	print("No such directory!")
	exit(1)

smean = loadtxt("%s/sall_set%d"%(directory,sett))
indices_L = loadtxt("%s/indicesL_set%d"%(directory,sett))

tm = loadtxt("%s/mall_set%d"%(directory,sett))[:,0]
mmean = loadtxt("%s/mall_set%d"%(directory,sett))[:,1:]

#-----------------------------plot data-------------------------------------------------------------#
  
fig, axs = plt.subplots(1, 1, figsize=(4,0.5))

cm = plt.get_cmap('gist_rainbow')
mycolors=[cm(1.*i/(p)) for i in range(p)]
#random.Random(1990).shuffle(mycolors)
mycolors=['#E6194B','#4363d8','#ffe119','#3cb44b','#f58231','#911eb4','#42d4f4','#f032e6','#bfef45','#469990', 'fuchsia', 'yellowgreen', 'magenta', 'orangered', 'crimson', 'aqua']

i=0
j=0
for mu in range(p):
	if mu in indices_L:
		axs.plot(tm, mmean[:,mu], color=mycolors[i], lw=0.5)
		i=i+1
	else:
		axs.plot(tm, mmean[:,mu], ls='--', lw=0.5)
		j=j+1

axs.set_ylabel('Overlap')
axs.set_xlabel('Time [sec]')
axs.set_xlim(0,50)	
axs.set_ylim(0,1.1)	
fig.savefig("m_ft_set%d_N%d_a%.2f_p%d_L%d_dJ%.4f_JEI%.4f.png"%(sett, N,a,p,L,dJ,J_EI), bbox_inches='tight')
#plt.show()