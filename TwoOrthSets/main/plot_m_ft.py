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
kappa=1.2

N=5000
a=0.01

p=600
U=0.5
Z=500
L=16
sett=0
pattset=0
T=500000

directory='output_N%d_a%.2f_p%d_L%d_dJ%.4f_kappa%.1f_set%d'%(N,a,p,L,dJ,kappa,pattset)

print(directory)

if os.path.isdir(directory) == False:
	print("No such directory!")
	exit(1)

smean = loadtxt("%s/sall_set%d"%(directory,sett))
indices_L = loadtxt("%s/indicesL_set%d"%(directory,sett))

tm = loadtxt("%s/mall_set%d"%(directory,sett))[:,0]
mmean = loadtxt("%s/mall_set%d"%(directory,sett))[:,1:]

#-----------------------------plot data-------------------------------------------------------------#
  
fig, axs = plt.subplots(1, 1, figsize=(6,1))

cm = plt.get_cmap('gist_rainbow')
mycolors=[cm(1.*i/(p)) for i in range(p)]
random.Random(1990).shuffle(mycolors)

for mu in range(p):
	if (mu in indices_L):
		axs.plot(tm, mmean[:,mu], color=mycolors[mu])
	else:
		axs.plot(tm, mmean[:,mu], color=mycolors[mu], ls='--')

axs.set_ylabel('Overlap')
axs.set_xlabel('Time [sec]')	
fig.savefig("m_ft_set%d_N%d_a%.2f_p%d_L%d_dJ%.4f_kappa%.1f.png"%(sett, N,a,p,L,dJ,kappa), bbox_inches='tight')
