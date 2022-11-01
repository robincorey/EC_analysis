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
import re

def get_data(res,pdb):
	print pdb
	for num in ['1','2','3','4','5']:
		print num
		file = open('Residue_distribution_CDL/z_%s_%s.pdb' % (pdb,num), "r")	
		for line in file:
			if re.search(res,line):
				parts = line.rstrip("\n").split(",")
				if float(parts[2]) > 0.1: 
					z.append(float(parts[1]))

def plot_data(z,res,i):
	n, bins, patches = plt.hist(z, bins=25,range=[-30,30], alpha=0)
	bin_mid = (bins[:-1] + bins[1:]) / 2
	plt.plot(bin_mid/10,n,color=colors[i], label=res) #,linestyle='dashed')
	limits.append(np.max(n))

colors = ['red','dodgerblue','darkorange','green','saddlebrown','lightblue']

pdbs = ['5OQT','1KQF','2IC8','5ZUG','4JR9','3O7P','1NEK','3ZE3','5MRW','5OC0','1ZCD','1PV6','3OB6','4IU8','1KF6','3QE7','1KPK','5SV0','1U77','1FFT','5AJI','4ZP0','5AZC','1Q16','1FX8','6AL2' , '2QFI','3DHW','1IWG','3K07','4GD3','5JWY','3B5D','4Q65','1PW4','2R6G','4DJI','2WSX','1L7V','1RC2']

fig = plt.figure()
ax = fig.add_subplot(111)

limits = []

for i in np.arange(1,20):
	try:
		res = sys.argv[i]
		z = []
        	for pdb in pdbs:
			get_data(res,pdb)
		plot_data(z, res,i-1)
	except IndexError:
		pass

plt.ylabel('lipid contacts', fontsize=25 )
plt.xlabel('z-axis position (nm)', fontsize=25 )
ax.plot([-2,-2],[0,2000],linewidth=15,color='gray',alpha=0.3)
ax.plot([2,2],[0,2000],linewidth=15,color='gray',alpha=0.3)
plt.xlim(-2.5,2.5)
plt.legend(fontsize=15,frameon=False)
plt.ylim(0,np.max(limits)+(0.1*np.max(limits)))
plt.savefig('Figs/Z_analysis_%s.cdl.highres.pdf' % sys.argv[1],bbox_inches='tight', dpi=1800,format='pdf')
