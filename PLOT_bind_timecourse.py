#!/sansom/s137/bioc1535/programs/anaconda2/bin/python

import scipy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys
from scipy import signal

plt.xlabel('time ($\mu$s)', fontsize=25 ) 
plt.ylabel('bound lipid count', fontsize=25)

# not yet updated

def plot_function (lipid, color):
	leaflets = [  'lower', 'upper' ]
	for shift, leaflet in enumerate(leaflets):
		import glob
		x = np.genfromtxt(fname='kinetics/upper/bind_5ZUG_1_CARD.xvg', usecols=(0),unpack=True)
		filename = glob.glob('kinetics/%s/bind_*_%s.xvg' % ( leaflet, lipid ))
		y = np.array([np.genfromtxt(f, usecols=(1),unpack=True) for f in filename])
		x_all = []
		y_all = []
		y_err = []
		for i in np.arange(2000):
			x_all.append(x[i]/1000)
			y_all.append(np.mean([lis[i] for lis in y]))
	#	y_smooth = sc.signal.savgol_filter(y_all, 11, 1, deriv=0, delta=1, axis=-1, mode='interp', cval=0.0)
		plt.plot(x_all, y_all, color=color, linewidth=1) #, alpha=shift/2)
	#	plt.plot(x_all, y_smooth, color=color, linewidth=1)
#	plt.savefig('%s_kinetics.highres.pdf' % (lipid), bbox_inches='tight', dpi=1800,format='pdf')

color = [ 'red', 'light' ]

fig = plt.figure()
ax = fig.add_subplot(111)
	
lipids = { "CARD":"royalblue" , "POPG":"red", "POPE":"orange" }

for lipid, color in lipids.items():
	plot_function(lipid, color)
	
plt.savefig('all_kinetics.highres.pdf', bbox_inches='tight', dpi=1800,format='pdf')
