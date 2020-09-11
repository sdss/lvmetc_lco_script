# Read ZP files in Francesco's format and write them in LCOETC format
# Mainly taking the maximum ZP in order overlap regions
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy
from astropy.io import ascii
from scipy.interpolate import interp1d

filein=sys.argv[1]
fileout=sys.argv[2]

data=ascii.read(filein, header_start=2, data_start=3)

order=data['ap_id'].data
lam=data['lambda'].data
zplam=data['zero_point'].data

uorder=np.unique(order)

lamout=np.linspace(lam.min(), lam.max(), len(lam))
zpaux=np.zeros([len(uorder), len(lamout)])


for i in range(len(uorder)):
    f=interp1d(lam[order == uorder[i]], zplam[order == uorder[i]], bounds_error=False, fill_value=0)
    zpaux[i,:]=f(lamout)


zpout=np.zeros(len(lamout))

for i in range(len(lamout)):
    zpout[i]=np.amax(zpaux[:,i])


ascii.write([lamout, zpout], fileout, format='commented_header', names=['lambda', 'zeropoint'], formats={'lambda':'%.3f', 'zeropoint':'%.2f'})

#plt.plot(lamout, zpaux[0,:])
plt.plot(lamout, zpout)
plt.show()


    
    
