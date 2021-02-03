import scipy as sc
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import re
import sys
from matplotlib.lines import Line2D


def get_total(pdb,leaflet):
	# how many lipids per pdb per leaflet
	all_lipids = []
	for lipid in lipids:
		lipid_total = [] # moved to here!
		for num in ['1','2','3','4','5']:
				ydata = np.loadtxt(fname='Leaflets/%s/%s/%s.%s.xvg' % ( leaflet,pdb,lipid,num ), usecols=(1), comments=['@','#'], unpack=True)
				lipid_total.append(np.mean(ydata))
		all_lipids.append(np.mean(lipid_total))
	return np.sum(all_lipids)					

def for_propensity(pdb,lipid,shift):
	# how many lipids are in leaflet
	top = '/sansom/s156a/bioc1535/Ecoli_patch/full_complement/chosen/%s/topol.top' % ( pdb)
	lipid_count = []
	for line in open(top, 'r'):
		if re.search(lipid, line):
			num = re.sub('[^0-9]', '', line)
			lipid_count.append(num)	
	return lipid_count[shift]

def plot_function(pdb, lipid, counter, shift, tot_lipids, leaflet,propensity):
		shift_value = float((shift - 0.5 )*0.35)
		x = [ counter+shift_value+(0.005*(np.random.randint(-5,5))) ]
		pdb_av = []
		for num in ['1', '2', '3','4','5']:
			ydata = np.loadtxt(fname='Leaflets/%s/%s/%s.%s.xvg' % ( leaflet,pdb,lipid,num ), usecols=(1), comments=['@','#'], unpack=True)
			pdb_av.append(np.mean((ydata/int(propensity)))) # different
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
        plt.text(x,np.max(alldata)+0.1,sig, horizontalalignment='center')
	return p2

leaflets = [ 'lower', 'upper' ]
lipids = [ 'CARD' , 'POPG', 'POPE' ]
pdbs = [ '5OQT'  ,'1KQF','2IC8','5ZUG','4JR9','3O7P','1NEK','3ZE3','5MRW','5OC0','1ZCD','1PV6','3OB6','4IU8','1KF6','3QE7','1KPK','5SV0','1U77','1FFT','5AJI','4ZP0','5AZC','1Q16','1FX8','6AL2', '2QFI','3DHW','1IWG','3K07','4GD3','5JWY','3B5D','4Q65','1PW4','2R6G','4DJI','2WSX','1L7V','1RC2']

fig = plt.figure()
ax = fig.add_subplot(111)
colors = [ 'royalblue', 'red', 'orange' ]
concs = [ 10, 23, 67 ]

for counter, lipid in enumerate(lipids):
	low = []
	up = []
	for shift, leaflet in enumerate(leaflets):
		for pdb in pdbs:
			tot_lipids = get_total(pdb,leaflet)
			# temp skipping # propensity = for_propensity(pdb,lipid,shift)
			propensity = concs[counter]
			plot_function(pdb, lipid, counter, shift, tot_lipids, leaflet,propensity)
	medianprops = dict(linewidth=2, color='black')
	bp = ax.boxplot(low, positions=[counter-0.175],showfliers=0,medianprops=medianprops)
	plt.setp(bp['boxes'], color='black',linewidth=2)
	bp2 = ax.boxplot(up, positions=[counter+0.175],showfliers=0,medianprops=medianprops)
	plt.setp(bp2['boxes'], color='black',linewidth=2)
	p2 = ttest(up,low,counter)
	print p2
	
plt.plot([-1,3.175],[1,1],color='gray', linestyle='dashed',linewidth=2,alpha=0.5)

plt.ylabel('lipid binding propensity', fontsize=25 )
plt.yticks(np.arange(0.5, 2.55, step=0.5),fontsize=20)
plt.xticks((-0.175,0.175,0.825,1.175,1.825,2.175), ('inner', 'outer', 'inner', 'outer', 'inner', 'outer'), rotation=45, fontsize=20)
plt.xlim([-0.5,2.5])
custom_lines = [Line2D([0], [0], color='royalblue', lw=2),
                Line2D([0], [0], color='red', lw=2),
                Line2D([0], [0], color='orange', lw=2)]
plt.legend(custom_lines, ['CDL','POPG','POPE'], fontsize=15,frameon=False, loc='upper right')
plt.savefig('Figs/Leaflets.propensity.png',bbox_inches='tight', dpi=900)
