import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import pylab
from scipy.optimize import fsolve

#see Mustill & Wyatt 2009 for a discussion of the equations used within...

debrisRingRadius = 110.0 #AU
stellarMass = 2.0 #solar masses
systemAge = 40e6 #yr

def eccen(m,a):
	return ((1.53e3*(debrisRingRadius/10.)**4.5*stellarMass**.5*(m)**(-1)*a**-3)/systemAge)**(2./3) ## eq 15 of Mustill and Wyatt (2009) solved for the eccentricity in (mass, semimajor axis) parameter space


#def e(a):
#	return debrisRingRadius**(1.5)*(3.8*(stellarMass)**(1./3)*(a)**(1./3)*(10)**(-2./3))**(-1.5) # eq 23 of Mustill & Wyatt (2009) solved for the eccentricity of the planet, e_{pl}

#def m(e,a):
#	return (1.53e3*(((1-e**2)**1.5)/e)*(a/10)**4.5*2**.5*a**-3)/(systemAge*1047.92612) ## eq 15 of M & W (2009) solved for the mass of the perturber required to stir the disk in a given time (in this case, the age of the system, ie 40Myr) as a function of the eccentricity, given above, and the semi major axis, a
# the factor of 1047.92612 comes from the fact that these eq are given in solar masses, we want jupiter masses to plot later


def make_cmap(colors, position=None, bit=False):
	import matplotlib as mpl
	bit_rgb = np.linspace(0,1,256)
	if position == None:
		position = np.linspace(0,1,len(colors))
	else:
		if len(position) != len(colors):
			sys.exit("position length must be the same as colors")
		elif position[0] != 0 or position[-1] != 1:
			sys.exit('position must start with 0 and end with 1')
	if bit:
		for i in range(len(colors)):
			colors[i] = (bit_rgb[colors[i][0]],					                             bit_rgb[colors[i][1]],					                             bit_rgb[colors[i][2]])
	cdict = {'red':[], 'green':[], 'blue':[]}
	for pos, color in zip(position, colors):
		cdict['red'].append((pos, color[0], color[0]))
		cdict['green'].append((pos, color[1], color[1]))
		cdict['blue'].append((pos, color[2], color[2]))
	cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	return cmap
colors = [(255,255,255), (95,245,255), (0,201,253), (0,56,233), (0,0,40)]
my_cmap = make_cmap(colors,bit=True)

rARR = np.arange(1,70.1,.3466) #radius array
mARR = np.arange(.1,10.1,0.05)*9.5458e-4 #mass array in solar masses, but we will label it as jupiter masses later


pixelMap = np.zeros([len(rARR),len(mARR)])
def imageGen(m,a):
	for i in range(len(rARR)):
		for j in range(len(mARR)):
			m,a = mARR[i],rARR[j]
			func = lambda e : (e**(2./3)/(1-e**2)) - eccen(m,a) #set up to solve for zeros
			#fsolve needs good hints depending on your location in parameter space in order to correctly solve this eq
			if a < 10 and a > 3 and m>1:
				pixelMap[i][j] = fsolve(func,.8)
			if a >3 and a < 10 and m <1:
                                pixelMap[i][j] = fsolve(func,.0)
			if a >0.5 and a <= 3 and m >1:
                                pixelMap[i][j] = fsolve(func,.95)
			if a >0.5 and a <= 3 and m <1:
                                pixelMap[i][j] = fsolve(func,.0)
			if a>=10 and a < 25 and m > 1:
				pixelMap[i][j] = fsolve(func,.5)  
			if a>=10 and a < 25 and m < 1:
				pixelMap[i][j] = fsolve(func,0)  
			if a>=25:
				pixelMap[i][j] = fsolve(func,0)  
			#if a>=50:
			#	pixelMap[i][j] = fsolve(func,0)  
imageGen(mARR,rARR)

#plot commands
plt.clf()
fig = plt.imshow(pixelMap,interpolation='nearest',cmap=my_cmap,origin='lower')
ax = pylab.gca()
formatter = matplotlib.ticker.FormatStrFormatter('%d')
ax.yaxis.set_major_formatter(formatter)
x=[1,10,20,30,40,50,60,70] 
y = [.1,2,4,6,8,10]
pylab.xticks([0,28,57,86,114,142,171,199],x,label='%d')
pylab.yticks([0,40,80,120,160,199],y,label='%d')
pylab.xlim(0,199)
pylab.ylim(0,199)
pylab.xlabel("$a_{pl}$",fontsize=20)
pylab.ylabel("$M_{pl} [M_{Jup}]$",fontsize=20)
manual_loc = [(14,60),(20,80),(25,100),(30,120),(37,140),(45,160),(60,180),(120,199)]
cset = plt.contour(pixelMap,[0.01,0.05,0.1,0.2,0.3,0.5,0.7,0.9],linewidths=2,fontsize=14,cmap=plt.cm.bone)
plt.clabel(cset,inline=True,fmt='%1.2f',fontsize=10,manual=manual_loc)
cbar = plt.colorbar(fig,ticks=[0.00,0.50,0.99]) # adding the colobar on the right
cbar.ax.set_yticklabels(['0.00','0.50','0.99'])
cbar.set_label("Eccentricity",rotation=90)
plt.savefig('49CET_planetStirring',dpi=400)
plt.show()
