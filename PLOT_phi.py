#!/sansom/s137/bioc1535/programs/anaconda2/bin/python

import scipy as sc
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys
import subprocess

def plot_function (): #num, res, lipid, color):
	card = []
	popg = []
	for i in np.arange(1,6):
		c,e,g = np.genfromtxt(fname='Phi_analysis/LIP.%s.txt' % i, delimiter=',', usecols=(1,3,5), unpack=True, skip_header=1)
		ax.scatter(np.arange(20),(c/np.sum(c))-(e/np.sum(e)),color='royalblue',s=20,alpha=0.5,linewidth=0)
		ax.scatter(np.arange(20),(g/np.sum(g))-(e/np.sum(e)),color='red',s=20,alpha=0.5, linewidth=0)
		card.append((c/np.sum(c))-(e/np.sum(e)))
		popg.append((g/np.sum(g))-(e/np.sum(e)))
	ax.errorbar(np.arange(20),np.mean(card,axis=0), marker='o', color='royalblue', linestyle='None',yerr=np.std(card,axis=0), linewidth=1, ecolor='royalblue', capsize=3)
	ax.errorbar(np.arange(20),np.mean(popg,axis=0), marker='o', color='red', linestyle='None',yerr=np.std(popg,axis=0), linewidth=1, ecolor='red', capsize=3) 

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot([0,21],[0,0], linewidth=2,linestyle='dashed', color='black')

plot_function() 

plt.ylabel('$\Delta(\Phi_{L}-\Phi_{POPE})$', fontsize=25 )
fig.set_size_inches(8, 5)
plt.xticks(np.arange(0, 20, step=1), ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'))
plt.savefig('Figs/phi.png', bbox_inches='tight', dpi=900)
