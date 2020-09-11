import sys
import numpy as np
import scipy
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy
from astropy.io import ascii
from astropy.table import Table


filein=sys.argv[1]
fileout=sys.argv[2]

# Convert radiance.dat file in units of micron vs photon/s/m2/micron/arcsec2
# to sky flux density file skyflambda.dat in units of A vs erg/s/cm2/A/arcsec2

#SYNTAX:
#> python rad2slam.py filein fileout


h=6.6260755e-27   # Planck's constant in [erg*s]
c=2.99792458e18   # Speed of light in [A/s]

# read radiance file
wavemu, rad = np.loadtxt(filein, unpack=1)

# convert units
wavea=wavemu*10.
slam=rad*h*c/wavea*1e-4*1e-4  # times h*nu times 1e-4 (per m2 to per cm2) times 1e-4 (per micron to per A)

# abmag=-2.5*np.log10(3.34e4*wavea**2*slam)+8.9

# write flux density file
np.savetxt(fileout, np.transpose([wavea, slam]), fmt=['%.3f', '%.5e'])
