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


def plot_function(res,pdb,component):
		counts = []
		for num in ['1','2','3','4','5']:
			file = open('Predicting_site_POPE/raw/%s_%s_%s.csv' % (pdb,num,component), "r")	
			for line in file:
				if re.search(res,line):
					parts = line.rstrip("\n").split(",")
					count = 0
    					iterparts = iter(parts)
					next(iterparts)
					for part in iterparts:
   						count = count + float(part)
					counts.append(count)
		total = np.mean(counts)
		return total
	        
components = [ 'NH3','PO4','GL2GL1GL3','CD' ]
residues =  ['GLY','ALA','ILE','LEU','PRO','VAL','PHE','TRP','TYR','ASP','GLU','ARG','HIS','LYS','SER','THR','CYS','MET','ASN','GLN']
pdbs = ['5OQT','1KQF','2IC8','5ZUG','4JR9','3O7P','1NEK','3ZE3','5MRW','5OC0','1ZCD','1PV6','3OB6','4IU8','1KF6','3QE7','1KPK','5SV0','1U77','1FFT','5AJI','4ZP0','5AZC','1Q16','1FX8','6AL2' , '2QFI','3DHW','1IWG','3K07','4GD3','5JWY','3B5D','4Q65','1PW4','2R6G','4DJI','2WSX','1L7V','1RC2']

fig = plt.figure(figsize=(12, 3))
ax = fig.add_subplot(111)

for component in components:
	data = [] #np.zeros((2,20))
	for x, res in enumerate(residues):
		plot = []
		for pdb in pdbs:
			total = plot_function(res,pdb,component)
			plot.append(total)
		data.append([res,np.nanmean(plot),sc.stats.sem(plot,nan_policy='omit')])
	sorted_data = sorted(data, key=lambda tup: tup[1], reverse=True)
	for x,i in enumerate(np.arange(20)):
		if sorted_data[i][0] in ('ARG', 'LYS'):
			col = 'dodgerblue'
		elif sorted_data[i][0] in ('TRP', 'PHE', 'TYR'):
			col = 'darkorange'
		elif sorted_data[i][0] in ('ASP', 'GLU'):
			col = 'red'
		elif sorted_data[i][0] in ('THR', 'SER','HIS','ASN','GLN'):
                        col = 'forestgreen'
		elif sorted_data[i][0] in ('GLY'):
                        col = 'purple'
		else:
			col = 'yellow'
		plt.bar(x,sorted_data[i][1],color=col,yerr=sorted_data[i][2],capsize=2,edgecolor='k',linewidth=1)
	plt.ylabel('contacts', fontsize=25 )
	plt.xticks(np.arange(20),(sorted_data[0][0],sorted_data[1][0],sorted_data[2][0],sorted_data[3][0],sorted_data[4][0],sorted_data[5][0],sorted_data[6][0],sorted_data[7][0],sorted_data[8][0],sorted_data[9][0],sorted_data[10][0],sorted_data[11][0],sorted_data[12][0],sorted_data[13][0],sorted_data[14][0],sorted_data[15][0],sorted_data[16][0],sorted_data[17][0],sorted_data[18][0],sorted_data[19][0]), rotation=45, ha='right',rotation_mode="anchor", fontsize=25)
	plt.ylim([0,sorted_data[0][1]*1.2])
	plt.yticks(fontsize=25)
	plt.savefig('Figs/Predicted_site_%s_allres_POPE.highres.pdf' % component, bbox_inches='tight', dpi=1800,format='pdf')
	plt.cla()
