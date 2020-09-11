import numpy as np
import astropy
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as plt

hdulist=fits.open('filter_curves.fits')

udata=hdulist[1].data
lu=udata['wavelength']
tu=udata['respt']

gdata=hdulist[2].data
lg=gdata['wavelength']
tg=gdata['respt']

rdata=hdulist[3].data
lr=rdata['wavelength']
tr=rdata['respt']

idata=hdulist[4].data
li=idata['wavelength']
ti=idata['respt']

zdata=hdulist[5].data
lz=zdata['wavelength']
tz=zdata['respt']


plt.plot(lu, tu)
plt.plot(lg, tg)
plt.plot(lr, tr)
plt.plot(li, ti)
plt.plot(lz, tz)
axes=plt.gca()
axes.set_xlim([3000, 10000])
axes.set_ylim([-0.1, 1])

plt.show()

ascii.write([lu, tu], 'sdss_u.dat', names=['wave', 'trans'], format='commented_header', formats={'wave':'%.2f', 'trans':'%.2f'}, overwrite=True)

ascii.write([lg, tg], 'sdss_g.dat', names=['wave', 'trans'], format='commented_header', formats={'wave':'%.2f', 'trans':'%.2f'}, overwrite=True)

ascii.write([lr, tr], 'sdss_r.dat', names=['wave', 'trans'], format='commented_header', formats={'wave':'%.2f', 'trans':'%.2f'}, overwrite=True)

ascii.write([li, ti], 'sdss_i.dat', names=['wave', 'trans'], format='commented_header', formats={'wave':'%.2f', 'trans':'%.2f'}, overwrite=True)

ascii.write([lz, tz], 'sdss_z.dat', names=['wave', 'trans'], format='commented_header', formats={'wave':'%.2f', 'trans':'%.2f'}, overwrite=True)




