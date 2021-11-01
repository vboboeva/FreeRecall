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
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}, size=12)
rc('text', usetex=True)
rc('axes', edgecolor='black', linewidth=0.5)
rc('legend', frameon=False)
rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}'] # \boldmath

#-----------------------------load data-------------------------------------------------------------#
p=300
alpha=0.032
beta=0.020
gamma=0.005
zeta=0.001
tau_pre=5
tau_post=5
U=0.5
Z=500
J_EI=0.25
L=10
sett=0
a=0.01
T=500000

directory='output_p%d_L%d_alpha%.4f_beta%.4f_gamma%.4f_zeta%.4f_taupre%.1f_kappa%.3f_pattset0'%(p,L,alpha,beta,gamma,zeta,tau_pre,J_EI)

print(directory)

if os.path.isdir(directory) == False:
	print("No such directory!")
	exit(1)

t = loadtxt("%s/input_set%d"%(directory,sett))[:,0]
stimulus_id = loadtxt("%s/input_set%d"%(directory,sett))[:,1]
stimulus_amp = loadtxt("%s/input_set%d"%(directory,sett))[:,2]
back_amp = loadtxt("%s/input_set%d"%(directory,sett))[:,3]


indices_L1 = loadtxt("%s/indicesL_set%d"%(directory,sett))[0,:]
indices_L2 = loadtxt("%s/indicesL_set%d"%(directory,sett))[1,:]
print(indices_L1)
print(indices_L2)

tm = loadtxt("%s/mall_set%d"%(directory,sett))[:,0]
mmean = loadtxt("%s/mall_set%d"%(directory,sett))[:,1:]
ts = loadtxt("%s/sall_set%d"%(directory,sett))[:,0]
smean = loadtxt("%s/sall_set%d"%(directory,sett))[:,1:]

#-----------------------------plot data-------------------------------------------------------------#
  

fig, axs = plt.subplots(2, 1, figsize=(9,1.5))

cm = plt.get_cmap('gist_rainbow')
mycolors=[cm(1.*i/(p)) for i in range(p)]
random.Random(1990).shuffle(mycolors)

for mu in range(p):
	if (mu in indices_L1):
		axs[0].plot(t[numpy.where(stimulus_id==mu)], stimulus_amp[numpy.where(stimulus_id==mu)], color=mycolors[mu])
	elif (mu in indices_L2):
		axs[0].plot(t[numpy.where(stimulus_id==mu)], stimulus_amp[numpy.where(stimulus_id==mu)], color=mycolors[mu])

#print(back_amp)
#axs[0].plot(t, back_amp, color='Red')
#axs[0].set_ylabel('$I$') 

for mu in range(p):
	if (mu in indices_L1):
		#axs[0].plot(ts, smean[:,mu], color='Gray')
		axs[1].plot(tm, mmean[:,mu], color=mycolors[mu], lw=0.5)
	elif (mu in indices_L2):
		#axs[0].plot(ts, smean[:,mu], color='Gray')
		axs[1].plot(tm, mmean[:,mu], color=mycolors[mu], ls='--', lw=0.5)		
	else:
		#axs[0].plot(ts, smean[:,mu], color='Gray')
		axs[1].plot(tm, mmean[:,mu], color='Gray', lw=0.5)#mycolors[mu], ls='--')

axs[0].set_xticklabels([])
axs[1].set_ylim(0,1.1)
axs[0].set_xlim(0,140)
axs[1].set_xlim(0,140)
axs[0].set_ylabel('Input') 
axs[1].set_ylabel('Overlap') 
axs[1].set_xlabel('Time [sec]')	
fig.savefig("m_ftp_set%d_p%d_L%d_alpha%.4f_beta%.4f_gamma%.4f_zeta%.4f_taupre%.1f_JEI%.3f.png"%(sett,p,L,alpha,beta,gamma,zeta,tau_pre,J_EI), bbox_inches='tight')