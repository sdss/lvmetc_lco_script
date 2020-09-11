import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy
from astropy.io import ascii
from scipy.interpolate import interp1d


data0=ascii.read('transmission_0.dat')

lam=10.*data0['col1']
klam=-2.5*np.log10(data0['col2'])

lam=lam[data0['col2'] != 0]
klam=klam[data0['col2'] != 0]


ascii.write([lam, klam], 'klam.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.5f'})


dataMcd=ascii.read('extinc_McD.dat')
dataCTIO=ascii.read('ctio_atmos.dat')
dataMK=ascii.read('maunakea_atmos.dat')


plt.plot(lam, klam, 'b', linewidth=2, label='Paranal')
plt.plot(dataMcd['col1']*10., dataMcd['col2'], 'r', linewidth=2, label='McDonald')
plt.plot(dataCTIO['col1'], dataCTIO['col2'], 'g', linewidth=2, label='CTIO')
plt.plot(dataMK['col1']*10., dataMK['col2'], 'm', linewidth=2, label='Mauna Kea')
plt.legend()
plt.show()






