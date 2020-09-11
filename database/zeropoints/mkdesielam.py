import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import astropy
from astropy.io import ascii
import scipy
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

data=ascii.read('desielam.txt')

x=data['col1'].data
y=data['col2'].data

xarr=np.linspace(3600, 9800, 1000)

sp=UnivariateSpline(x,y)

sp.set_smoothing_factor(5e-3)

ascii.write([xarr, sp(xarr)], 'LVM_LCO25_ELAM.dat', names=['lambda', 'elam'], format='commented_header', overwrite=True, formats={'lambda':'%3f', 'elam':'%2f'})

ascii.write([xarr, sp(xarr)], 'LVM_APO25_ELAM.dat', names=['lambda', 'elam'], format='commented_header', overwrite=True, formats={'lambda':'%3f', 'elam':'%2f'})

ascii.write([xarr, sp(xarr)], 'LVM_LVM160_ELAM.dat', names=['lambda', 'elam'], format='commented_header', overwrite=True, formats={'lambda':'%3f', 'elam':'%2f'})

ascii.write([xarr, sp(xarr)], 'LVM_NMSU10_ELAM.dat', names=['lambda', 'elam'], format='commented_header', overwrite=True, formats={'lambda':'%3f', 'elam':'%2f'})

plt.plot(x,y)

plt.plot(xarr,sp(xarr))

plt.show()

