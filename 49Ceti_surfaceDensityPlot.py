from random import randint
from scipy.stats import mstats
from matplotlib.ticker import FormatStrFormatter, LogLocator
import numpy as np
from pyfits import getdata
import matplotlib.pyplot as plt
from disk_SPL import Disk as SPL
from disk_SPLWR import Disk as SPLWR
from disk_DPL import Disk as DPL

'''
SPL = single power law
SPLWR = single power law with ring
DPL = double power law
'''

num_walkers = 100 # number of walkers to grab surface density profiles from

# final MCMC runs of the various models
whatbywhat = '60x900'
SPLchainDat = getdata('MCMCRUNS/SPL/'+whatbywhat+'/'+whatbywhat+'.chain.fits')
SPLchiDat = getdata('MCMCRUNS/SPL/'+whatbywhat+'/'+whatbywhat+'.chi.fits')

whatbywhat = '100x1200'
SPLWRchainDat = getdata('MCMCRUNS/SPLWR/'+whatbywhat+'/'+whatbywhat+'.chain.fits')
SPLWRchiDat = getdata('MCMCRUNS/SPLWR/'+whatbywhat+'/'+whatbywhat+'.chi.fits')

whatbywhat = '40x800'
DPLchainDat = getdata('MCMCRUNS/DPL/'+whatbywhat+'/'+whatbywhat+'.chain.fits')
DPLchiDat = getdata('MCMCRUNS/DPL/'+whatbywhat+'/'+whatbywhat+'.chi.fits')

def findBestFit(diskChain,diskChi):
	bestFit = np.ravel(np.array(np.where(diskChi < np.min(diskChi)+.000001)))
	bestDisk = diskChain[bestFit[0],bestFit[len(bestFit)/2]]
	return bestDisk

SPL_bestFit = findBestFit(SPLchainDat,SPLchiDat) 
SPLWR_bestFit = findBestFit(SPLWRchainDat,SPLWRchiDat)
DPL_bestFit = findBestFit(DPLchainDat,DPLchiDat)

#SPL_best_disk = SPL(SPL_bestFit[0],SPL_bestFit[1],SPL_bestFit[2],10**SPL_bestFit[3],10**SPL_bestFit[4],10**SPL_bestFit[5],SPL_bestFit[6],SPL_bestFit[7])
#SPLWR_best_disk = SPLWR(SPLWR_bestFit[0],SPLWR_bestFit[1],SPLWR_bestFit[2],10**SPLWR_bestFit[3],10**SPLWR_bestFit[4],10**SPLWR_bestFit[5],10**SPLWR_bestFit[6],SPLWR_bestFit[7],SPLWR_bestFit[8],SPLWR_bestFit[9])
#DPL_best_disk = DPL(DPL_bestFit[0],DPL_bestFit[1],DPL_bestFit[2],DPL_bestFit[3],10**DPL_bestFit[4],10**DPL_bestFit[5],10**DPL_bestFit[6],DPL_bestFit[7],DPL_bestFit[8],DPL_bestFit[9])

def findRandomFit(diskChain,diskType): # pulls random model parameters corresponding to somewhere on the MCMC chain after the burn in phase
	randomDisk = diskChain[randint(0,len(diskChain[:,0,0]))-1,randint(300,len(diskChain[0,:,0]))-1,:] ## all runs are burned in by 300 steps
	return randomDisk

def makeRandomModels(diskChain,diskType):
	rmp = findRandomFit(diskChain,diskType)
	# the models and what parameters they call are set up slightly differently. accounting for that here. 
	if diskType == 'SPL':
		if rmp is not None:
			randomModelDisk = SPL(rmp[0],rmp[1],rmp[2],10**rmp[3],10**rmp[4],10**rmp[5],rmp[6],rmp[7])
		else: return None
	if diskType == 'SPLWR':
		if rmp is not None:
			randomModelDisk = SPLWR(rmp[0],rmp[1],rmp[2],10**rmp[3],10**rmp[4],10**rmp[5],10**rmp[6],rmp[7],rmp[8],rmp[9])  
		else: return None
	if diskType == 'DPL':
		if rmp is not None:
			randomModelDisk = DPL(rmp[0],rmp[1],rmp[2],rmp[3],10**rmp[4],10**rmp[5],10**rmp[6],rmp[7],rmp[8],rmp[9])    
		else: return None
	return randomModelDisk


def findSurfaceDensity(diskType,diskModel): #SPL_bestFit, SPL_best_disk
	radius = np.arange(0,350,.1)
	GDarr = np.zeros(len(radius))
	# calculate the surface density distribution, taking the sum of all components of the disk. the SPL and DPL models only have two parts (inner and outer disk) whereas the SPLWR model has these two parts plus the additional ring. 
	if diskType == 'SPL' or diskType == 'DPL':
		for i in range(len(radius)):
			GDarr[i] = (diskModel.calculateGrainDistribution(radius[i]*1.496e11)*diskModel.grainMass+diskModel.calculateGrainDistribution3(radius[i]*1.496e11)*diskModel.grainMass3)
	else:
		for i in range(len(radius)):
			GDarr[i] = (diskModel.calculateGrainDistribution(radius[i]*1.496e11)+diskModel.calculateGrainDistribution2(radius[i]*1.496e11))*diskModel.grainMass+diskModel.calculateGrainDistribution3(radius[i]*1.496e11)*diskModel.grainMass3
	return GDarr

# initialize empty walkers of size num_walkers by the disk width... here we are doing .1 to 350 AU at a resolution of 0.1AU, so 3500 steps
SPL_distr = np.array([[[0.] for x in range(num_walkers)] for x in range(3500)])
SPLWR_distr = np.array([[[0.] for x in range(num_walkers)] for x in range(3500)])
DPL_distr = np.array([[[0.] for x in range(num_walkers)] for x in range(3500)])

def plotSwarm():
	for i in range(num_walkers):
		print i # keep track of how long this has been running
		#SPLWR random model
		SPLWR_randomModel = makeRandomModels(SPLWRchainDat,'SPLWR')
		if SPLWR_randomModel is not None:
			SPLWR_xSurf = findSurfaceDensity('SPLWR',SPLWR_randomModel)
			for j in range(len(SPLWR_xSurf)):
				SPLWR_distr[j][i] = SPLWR_xSurf[j] # save the surface density every 0.1AU for each walker
			
		#DPL random model
		DPL_randomModel = makeRandomModels(DPLchainDat,'DPL')
		if DPL_randomModel is not None:
			DPL_xSurf = findSurfaceDensity('DPL',DPL_randomModel)
			for j in range(len(DPL_xSurf)):
				
				DPL_distr[j][i] = DPL_xSurf[j]

		#SPL random model
		SPL_randomModel = makeRandomModels(SPLchainDat,'SPL')
		if SPL_randomModel is not None:
			SPL_xSurf = findSurfaceDensity('SPL',SPL_randomModel)
			for j in range(len(SPL_xSurf)):
				SPL_distr[j][i] = SPL_xSurf[j]

plotSwarm()

def distribution(distr):
	lower_lim = np.zeros(3500)
	median = np.zeros(3500)
	upper_lim = np.zeros(3500)
	for r in range(len(distr)):
		radii_temp = np.ravel(distr[r,:])
		quantiles = mstats.mquantiles(radii_temp,prob=[.16,.5,.84],axis=None)
		# ^ for a given radius, what is the distribution of surface densities for a given model? save .16, .5, and .84 for -1sigma, median, and +1sigma
		lower_lim[r] = quantiles[0]
		median[r] = quantiles[1]
		upper_lim[r] = quantiles[2]
	return lower_lim,median,upper_lim

#unpack distributions into variables for easy plotting
DPL_lower_lim,DPL_median,DPL_upper_lim = distribution(DPL_distr)
SPLWR_lower_lim,SPLWR_median,SPLWR_upper_lim = distribution(SPLWR_distr)
SPL_lower_lim,SPL_median,SPL_upper_lim = distribution(SPL_distr)


def plotSurfaceDensity():
	plt.clf()
	axes = plt.subplot(111)
	x = np.arange(0,350,.1)
	# each model has a dotted line for where we do not actually resolve the surface density profile (<40AU for DPL, <73AU for SPL, <60AU for SPLWR) and is solid thereafter. the fill denotes the region within 1sigma of the median value for each model. 

	plt.loglog(np.arange(np.min(np.where(DPL_lower_lim!=0))*.1,39.95,.1),1e4*DPL_median[np.where((DPL_lower_lim!=0) & (x<40))],color='r',linestyle=':',linewidth=3)
	plt.loglog(np.arange(40,np.max(np.where(DPL_lower_lim!=0))*.1+.05,.1),1e4*DPL_median[np.where((DPL_lower_lim!=0) & (x>39.95))],color='r',linewidth=3,label='Double Power Law')
	plt.fill_between(np.arange(np.min(np.where(DPL_lower_lim!=0))/10.,np.max(np.where(DPL_lower_lim!=0))/10.+.05,.1),1e4*DPL_lower_lim[np.where(DPL_lower_lim!=0)],1e4*DPL_upper_lim[np.where(DPL_lower_lim!=0)],color='r',alpha='.3')
		
	plt.loglog(np.arange(np.min(np.where(SPL_lower_lim!=0))/10.,72.95,.1),1e4*SPL_median[np.where((SPL_lower_lim!=0) & (x < 73))],color='b',linestyle=':',linewidth=3)
	plt.loglog(np.arange(73,np.max(np.where(SPL_lower_lim!=0))/10.+.05,.1),1e4*SPL_median[np.where((SPL_lower_lim!=0) & (x > 72.95))],color='b',linewidth=3,label='Single Power Law')
	plt.fill_between(np.arange(np.min(np.where(SPL_lower_lim!=0))/10.,np.max(np.where(SPL_lower_lim!=0))/10.+.05,.1),1e4*SPL_lower_lim[np.where(SPL_lower_lim!=0)],1e4*SPL_upper_lim[np.where(SPL_lower_lim!=0)],color='b',alpha='.3')

	plt.loglog(np.arange(np.min(np.where(SPLWR_lower_lim!=0))/10.,59.95,.1),1e4*SPLWR_median[np.where((SPLWR_lower_lim!=0) & (x<60))],color='g',linewidth=3,linestyle=':')	
	plt.loglog(np.arange(60,np.max(np.where(SPLWR_lower_lim!=0))/10.+0.05,.1),1e4*SPLWR_median[np.where((SPLWR_lower_lim!=0) & (x > 59.95))],color='g',linewidth=3,label='Single Power Law With Ring')
	plt.fill_between(np.arange(np.min(np.where(SPLWR_lower_lim!=0))/10.,np.max(np.where(SPLWR_lower_lim!=0))/10.+.05,.1),1e4*SPLWR_lower_lim[np.where(SPLWR_lower_lim!=0)],1e4*SPLWR_upper_lim[np.where(SPLWR_lower_lim!=0)],color='g',alpha='.3')

	# old, without the dotted/undotted bit
	#plt.loglog(np.arange(np.min(np.where(DPL_lower_lim!=0))/10.,np.max(np.where(DPL_lower_lim!=0))/10.+.05,.1),DPL_median[np.where(DPL_lower_lim!=0)])	
	#plt.fill_between(np.arange(np.min(np.where(DPL_lower_lim!=0))/10.,np.max(np.where(DPL_lower_lim!=0))/10.+.05,.1),DPL_lower_lim[np.where(DPL_lower_lim!=0)],DPL_upper_lim[np.where(DPL_lower_lim!=0)])
	#plt.loglog(np.arange(np.min(np.where(SPL_lower_lim!=0))/10.,np.max(np.where(SPL_lower_lim!=0))/10.+.05,.1),SPL_median[np.where(SPL_lower_lim!=0)])	
	#plt.fill_between(np.arange(np.min(np.where(SPL_lower_lim!=0))/10.,np.max(np.where(SPL_lower_lim!=0))/10.+.05,.1),SPL_lower_lim[np.where(SPL_lower_lim!=0)],SPL_upper_lim[np.where(SPL_lower_lim!=0)])
	#plt.loglog(np.arange(np.min(np.where(SPLWR_lower_lim!=0))/10.,np.max(np.where(SPLWR_lower_lim!=0))/10.+0.05,.1),SPLWR_median[np.where(SPLWR_lower_lim!=0)])	
	#plt.fill_between(np.arange(np.min(np.where(SPLWR_lower_lim!=0))/10.,np.max(np.where(SPLWR_lower_lim!=0))/10.+.05,.1),SPLWR_lower_lim[np.where(SPLWR_lower_lim!=0)],SPLWR_upper_lim[np.where(SPLWR_lower_lim!=0)])

	plt.fill_between(np.arange(20,30.1,.1),1e-3,1e2,color='k',alpha='.1') # shaded box denoting the resolution of our data
	axes.set_xticks(np.array([20,30,40,50,60,70,80,90,100,200,300]))	
	axes.set_xticklabels(['20','30','40','50','','70','','','100','200','300'])
	for tick in axes.xaxis.get_major_ticks():
		tick.label.set_fontsize(16) 
 	for tick in axes.xaxis.get_minor_ticks():
		tick.label.set_fontsize(16) 
 	for tick in axes.yaxis.get_major_ticks():
		tick.label.set_fontsize(16) 
 	for tick in axes.yaxis.get_minor_ticks():
		tick.label.set_fontsize(16) 
	plt.xlabel('Radius (AU)',fontsize=20)
	plt.ylabel('Surface Density (g/cm$^{2}$)', fontsize=20)
	plt.xlim(20,320)
	plt.ylim(3e-3,10)
	#plt.savefig('/home/jliemansifry/Desktop/49CET_SDD.png',dpi=400)
	plt.show()

plotSurfaceDensity()
