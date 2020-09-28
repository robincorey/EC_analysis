import scipy as sc
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys
from matplotlib.lines import Line2D


def get_total(pdb,leaflet):
	all_lipids = []
	for lipid in lipids:
		lipid_total = [] # moved to here!
		for num in ['1','2','3','4','5']:
				ydata = np.loadtxt(fname='Leaflets/%s/%s/%s.%s.xvg' % ( leaflet,pdb,lipid,num ), usecols=(1), comments=['@','#'], unpack=True)
				lipid_total.append(np.mean(ydata))
		all_lipids.append(np.mean(lipid_total))
	return np.sum(all_lipids)					
	
def plot_function(pdb, lipid, counter, shift, tot_lipids, leaflet):
		shift_value = float((shift - 0.5 )*0.35)
		x = [ counter+shift_value+(0.005*(np.random.randint(-5,5))) ]
		pdb_av = []
		for num in ['1', '2', '3','4','5']:
			ydata = np.loadtxt(fname='Leaflets/%s/%s/%s.%s.xvg' % ( leaflet,pdb,lipid,num ), usecols=(1), comments=['@','#'], unpack=True)
			pdb_av.append(np.mean(ydata))
		pdb_tot = np.mean(pdb_av)
		if leaflet is 'lower':
			low.append((pdb_tot/tot_lipids)*100)
		elif leaflet is 'upper':
			up.append((pdb_tot/tot_lipids)*100)
		ax.scatter(x, (pdb_tot/tot_lipids)*100, facecolors='none', s=25, edgecolors=colors[counter], linewidths=1, alpha=0.5) #0.2)
	        
def ttest(up, low, x): 
	p2 = 1
        t2, p2 = stats.ttest_ind(up,low)
        if p2 < 100:
                sig = 'ns'
        if p2 < 0.05:
                sig = '*'
        if p2 < 0.01:
                sig = '**'
        if p2 < 0.001:
                sig = '***'
        alldata = np.concatenate((up, low))
        plt.text(x,np.max(alldata)+1,sig, horizontalalignment='center')
	return p2

leaflets = [ 'lower', 'upper' ]
lipids = [ 'POPE' , 'POPG', 'CARD' ]
pdbs = [ '5OQT'  ,'1KQF','2IC8','5ZUG','4JR9','3O7P','1NEK','3ZE3','5MRW','5OC0','1ZCD','1PV6','3OB6','4IU8','1KF6','3QE7','1KPK','5SV0','1U77','1FFT','5AJI','4ZP0','5AZC','1Q16','1FX8','6AL2', '2QFI','3DHW','1IWG','3K07','4GD3','5JWY','3B5D','4Q65','1PW4','2R6G','4DJI','2WSX','1L7V','1RC2']

fig = plt.figure()
ax = fig.add_subplot(111)
colors = [ 'orange', 'red', 'royalblue' ]

for counter, lipid in enumerate(lipids):
	low = []
	up = []
	for shift, leaflet in enumerate(leaflets):
		for pdb in pdbs:
			tot_lipids = get_total(pdb,leaflet)
			plot_function(pdb, lipid, counter, shift, tot_lipids, leaflet)
	medianprops = dict(linewidth=2, color='black')
	bp = ax.boxplot(low, positions=[counter-0.175],showfliers=0,medianprops=medianprops)
	plt.setp(bp['boxes'], color='black',linewidth=2)
	bp2 = ax.boxplot(up, positions=[counter+0.175],showfliers=0,medianprops=medianprops)
	plt.setp(bp2['boxes'], color='black',linewidth=2)
	p2 = ttest(up,low,counter)
	print p2
	
plt.plot([-1,3.175],[67,67],color='orange', linestyle='dashed',linewidth=1,alpha=0.5)
plt.plot([-1,3.175],[23,23],color='red', linestyle='dashed',linewidth=1,alpha=0.5)
plt.plot([-1,3.175],[10,10],color='royalblue',linestyle='dashed', linewidth=1,alpha=0.5)

plt.ylabel('lipids bound (%)', fontsize=25 )
plt.yticks(np.arange(0, 71, step=10),fontsize=20)
plt.xticks((-0.175,0.175,0.825,1.175,1.825,2.175), ('inner', 'outer', 'inner', 'outer', 'inner', 'outer'), rotation=45, fontsize=20)
plt.xlim([-0.5,2.5])
custom_lines = [Line2D([0], [0], color='orange', lw=2),
                Line2D([0], [0], color='red', lw=2),
                Line2D([0], [0], color='royalblue', lw=2)]
plt.legend(custom_lines, ['POPE','POPG','CDL'], fontsize=15,frameon=False, loc='upper right')
plt.savefig('Figs/Leaflets.png',bbox_inches='tight', dpi=900)
