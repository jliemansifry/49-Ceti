## Notes related to the 49 Ceti project...

## Background: 
49 Ceti is a gas-rich debris disk that was recently resolved by ALMA at unprecedented angular resolution and 
sensitivity. My senior thesis involved modeling these data in order to better understand the processes that shape planet 
formation. We find that the dust disk of 49 Ceti is unlike anything observed before. It has two distinct regions, with an inner disk characterized by small , 0.1 micron grains, whereas the outer belt is characterized by larger, 2.0 micron grains. The disk
increases in surface density to a peak around 110AU before falling off thereafter. We think this region of higher density in
the dust corresponds to a similar increase in the gas surface density, but models have thus far been inconclusive.

## Data presentation: 
49Ceti_dataPlots.py takes 49Ceti_.75arcsectaper.fits and 49Ceti_robust2.fits as inputs and creates 
49Ceti_dataPlots.png as output. The color maps used are custom in order to appear continuous in color and greyscale. 
It could be easily edited to take different data images if desired.

## What is responsible for the peak in surface density? 
49Ceti_plotStirring.py generates a plot (49 Ceti_plotStirring.png)
of the eccentricity needed for a hypothetical planet in semi-major axis/mass space in order to "stir" the planetesimal belt at 
110AU. Debris rings need something to jostle them to continually replenish the small dust grains we observe. 49 Ceti's disk
is not massive enough to have "self-stirred" this ring in the age of the system, and the only other option is planetary
stirring. Using formulas from Mustill & Wyatt (2009), we are able to get a sense of the combination of the
eccentricity/mass/semimajor necessary to stir the ring. One such example is that a Jupiter-like planet (eccentricity ~ 0.05, 
mass ~ 1M_Jup) would need an orbit with a semimajor axis of 40AU to stir the ring. See the plot for all possible combinations
of parameters. 
