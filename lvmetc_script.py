# Local Volume Mapper (LVM) Exposure Time Calculator
# adapted from the Las Campanas Observatory ETC
#
# 2017-08-08: This is a BETA Version written by Guillermo A. Blanc. It has not been validated. Use at your own risk.
# 2020-09-09: Ported from a cgi web interface to a stand alone python script by Kathryn Kreckel.



import os,sys

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import astropy
from astropy.io import ascii
from astropy.io import fits
import astropy.convolution as conv
import scipy
from scipy.interpolate import interp1d
from datetime import datetime

# ==============================================
# ==============================================
# FUNCTION DEFINITIONS


def flam2epp(lam, flam, ddisp):
    # Convert flux density [erg/s/cm2/A] to photons per pixel [photons/s/cm2/pixel]
    # ddisp is the pixel scale in [A/pixel]
    h=6.6260755e-27   # Planck's constant in [erg*s]
    c=2.99792458e18   # Speed of light in [A/s]
    return flam*lam*ddisp/(h*c)

def zp2elam(zp, lam, neobj, ddisp, texp, atel, klam, amass): # System Efficiency from Zero Point
    # Claculates the system efficiency "elam" from a Zero Point "zp" (AB mag) corresponding to "neobj" 
    # electrons per pixel of width "ddisp" collected in "texp" at airmass "amass"
    # The ZP is typically given for 1 e/s/pixel at amass=1 and sometimes 1 ADU/s/pixel at amass=1
    # !!!!!!!! In the latter case neobj should be gain*naduobj ????????? 
    h=6.6260755e-27   # Planck's constant in [erg*s]
    return h*lam*neobj/(ddisp*texp*atel)*10**(0.4*(zp+48.6+klam*amass))

def elam2zp(elam, lam, neobj, ddisp, texp, atel, klam, amass): # Zero Point from System Efficiency
    # inverse function of zp2elam
    h=6.6260755e-27   # Planck's constant in [erg*s]
    return 2.5*np.log10(elam*ddisp*texp*atel/(h*lam*neobj))-48.6-klam*amass

def mklam(lrange, ddisp):
    # Returns wavelength array sampled at the instrument pixel dispersion
    return np.arange(lrange[0], lrange[1]+ddisp, ddisp)
    
def mkklam(lam, dlsf, dfiblam, klamfile):
    # Returns atmospheric extinction coefficient sampled at wavelengths in "lam"
    # The klam curve is convolved to the spectrograph LSF and the fiber width before sampling it
    #data=ascii.read('/Users/gblancm/work/LCO/lcoetc/database/sky/klam.dat')
    data=ascii.read(klamfile)    
    # trim to lrange only
    lam0=data['col1'][(data['col1'] > np.amin(lam)) * (data['col1'] < np.amax(lam))]
    klam0=data['col2'][(data['col1'] > np.amin(lam)) * (data['col1'] < np.amax(lam))]
    # interpolate to regular grid
    ddisp0=np.median(lam0[1:len(lam0)-1]-lam0[0:len(lam0)-2])
    lam1=np.arange(np.amin(lam0), np.amax(lam0), ddisp0)
    f=interp1d(lam0, klam0, fill_value='extrapolate')
    klam1=f(lam1)
    # convolve with LSF (Gaussian of FWHM=dlsf)
    gauss=conv.Gaussian1DKernel(stddev=dlsf/2.355/ddisp0, x_size=round_up_to_odd(10*dlsf/2.355/ddisp0))
    klam2=conv.convolve(klam1, gauss.array, boundary='extend')
    # convolve with fiber profile (Gaussian of FWHM=dfiberlam)
    gauss=conv.Gaussian1DKernel(stddev=dfiblam/2.355/ddisp0, x_size=round_up_to_odd(10*dfiblam/2.355/ddisp0))
    klam3=conv.convolve(klam2, gauss.array, boundary='extend')
    # interpolate to lam
    f=interp1d(lam1, klam3, fill_value='extrapolate')
    return f(lam)

def mkslam(lam, dlsf, dfiblam, sfile):
    # Returns sky surface brightness spectrum for a given moon phase, LSF FWHM, and slit width.
    # The result is convolved with LSF+slit and sampled at wavelengths in lam
    data=ascii.read(sfile)
    # trim to lrange only
    lam0=data['col1'][(data['col1'] > np.amin(lam)) * (data['col1'] < np.amax(lam))]
    slam0=data['col2'][(data['col1'] > np.amin(lam)) * (data['col1'] < np.amax(lam))]
    # interpolate to regular grid
    ddisp0=np.median(lam0[1:len(lam0)-1]-lam0[0:len(lam0)-2])
    lam1=np.arange(np.amin(lam0), np.amax(lam0), ddisp0)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    # convolve with LSF (Gaussian of FWHM=dlsf)
    gauss=conv.Gaussian1DKernel(stddev=dlsf/2.355/ddisp0, x_size=round_up_to_odd(10*dlsf/2.355/ddisp0))
    slam2=conv.convolve(slam1, gauss.array, boundary='extend')
    # convolve with fiber profile (Gaussian with FWHM=dfiblam)
    gauss=conv.Gaussian1DKernel(stddev=dfiblam/2.355/ddisp0, x_size=round_up_to_odd(10*dfiblam/2.355/ddisp0))
    slam3=conv.convolve(slam2, gauss.array, boundary='extend')
    logf.write(str(np.mean(slam1))+'\n')
    logf.write(str(np.mean(slam2))+'\n')
    logf.write(str(np.mean(slam3))+'\n')
    # interpolate to lam
    f=interp1d(lam1, slam3, fill_value='extrapolate')
    return f(lam)

def mkolam(lam, dlsf, dfiblam, ofile, filfile, abmag):
    # Returns object "ofile" spectrum for a given LSF FWHM, and slit width.
    # The result is sampled in wavelength at the instrument pixel dispersion
    # Result does not include slitlosses. Normalized to ABmag in filter in "filfile"
    c=2.99792458e18   # Speed of light in [A/s]
    data=ascii.read(ofile, names=['col1', 'col2'])
    lam0=data['col1']
    olam0=data['col2']
    # interpolate to regular grid (finer between original file grid and lam)
    ddisp0=np.amin([np.median(lam0[1:len(lam0)-1]-lam0[0:len(lam0)-2]), np.median(lam[1:len(lam)-1]-lam[0:len(lam)-2])])
    lam1=np.arange(np.amin(lam0), np.amax(lam0), ddisp0)
    f=interp1d(lam0, olam0, fill_value='extrapolate')
    olam1=f(lam1)
    # convolve with LSF (Gaussian of FWHM=dlsf)
    gauss=conv.Gaussian1DKernel(stddev=dlsf/2.355/ddisp0, x_size=round_up_to_odd(10*dlsf/2.355/ddisp0))
    olam2=conv.convolve(olam1, gauss.array, boundary='extend')
    # convolve with fiber profile (Gaussian with FWHM=dfiblam)
    gauss=conv.Gaussian1DKernel(stddev=dfiblam/2.355/ddisp0, x_size=round_up_to_odd(10*dfiblam/2.355/ddisp0))
    olam3=conv.convolve(olam2, gauss.array, boundary='extend')
    # read filter curve and interpolate to object lambda grid
    datafil=ascii.read(filfile, names=['col1', 'col2'])
    lamfil0=datafil['col1']
    tfil0=datafil['col2']
    f=interp1d(lamfil0, tfil0, fill_value='extrapolate')
    tfil1=f(lam1)
    tfil1[lam1<=lamfil0.min()]=0 
    tfil1[lam1>=lamfil0.max()]=0
    lameff=np.trapz(lam1*tfil1*lam1, lam1)/np.trapz(tfil1*lam1, lam1)
    # normalize to ABmag in filter "filfile"
    monoflam0=np.trapz(olam3*tfil1*lam1, lam1)/np.trapz(tfil1*lam1, lam1)
    monoflam1=10.0**(-0.4*(abmag+48.6))*c/lameff**2
    olam4=olam3*monoflam1/monoflam0
    logf.write(str(monoflam0)+'\n')
    logf.write(str(monoflam1)+'\n')
    logf.write(str(np.mean(olam1))+'\n')
    logf.write(str(np.mean(olam2))+'\n')
    logf.write(str(np.mean(olam3))+'\n')
    logf.write(str(np.mean(olam4))+'\n')
    # interpolate to lam
    f=interp1d(lam1, olam4, fill_value='extrapolate')
    return f(lam)

def mklinelam(lam, dlsf, dfiblam, linelam, lineflux, linefwhm):
    sigma=linefwhm/2.355
    linelam1=lineflux/(np.sqrt(2)*sigma)*np.exp(-1*(lam-linelam)**2/(2*sigma**2))
    # convolve with LSF (Gaussian of FWHM=dlsf)
    ddisp0=np.median(lam[1:len(lam)-1]-lam[0:len(lam)-2])
    gauss=conv.Gaussian1DKernel(stddev=dlsf/2.355/ddisp0, x_size=round_up_to_odd(10*dlsf/2.355/ddisp0))
    linelam2=conv.convolve(linelam1, gauss.array, boundary='extend')
    # convolve with fiber profile (Gaussian with FWHM=dfiblam)
    gauss=conv.Gaussian1DKernel(stddev=dfiblam/2.355/ddisp0, x_size=round_up_to_odd(10*dfiblam/2.355/ddisp0))
    linelam3=conv.convolve(linelam2, gauss.array, boundary='extend')
    return linelam3


def mkelam(lam, elamfile):
    # Returns efficiency curve from elam file (must be fractional transmission)
    data=ascii.read(elamfile, names=['col1', 'col2'])
    lam0=data['col1']
    elam0=data['col2']
    f=interp1d(lam0, elam0, fill_value='extrapolate')
    elam1=f(lam)
    return elam1

def mkneobj(lam, olam, ddisp, texp, elam, atel, klam, amass): 
    # Number of object electrons per spectral pixel
    # Returns the total number of electrons in a spectral pixel of width ddisp. 
    # Assumes all electrons fall in one spatial pixel so the output must then be redistribute across the spatial PSF.
    return flam2epp(lam, olam, ddisp)*texp*elam*atel*10**(-0.4*(klam*amass))

def mkneobj2d(neobj, dfib, nypix): 
    # Return a 2D spectrum of object (len(lam) x nypix)
    # nypix gets ceil to the next odd integer
    nypix=int(round_up_to_odd(nypix))
    ne2d=np.zeros((len(neobj),nypix))
    nefrac=np.zeros(nypix)
    for i in range(nypix):
        j=i-int(np.floor(nypix/2))
        ymin=(j-0.5)/np.sqrt(2)/(dfib/2.355)
        ymax=(j+0.5)/np.sqrt(2)/(dfib/2.355)
        nefrac[i]=(scipy.special.erf(ymax)-scipy.special.erf(ymin))/2. 
               
    for i in range(len(neobj)):
        ne2d[i,:]=neobj[i]*nefrac
        
    return ne2d
    
def mknesky(lam, slam, ddisp, texp, elam, atel):
    # Number of electrons per spectral pixel
    # Returns the number of electrons in one spatial pixel ddisp given a slit width dslit and a sky 
    # surface brightness spectrum slam.
    return flam2epp(lam, slam, ddisp)*texp*elam*atel

def mknesky2d(nesky, dfib, nypix):
    # Return a 2D spectrum of the sky (len(lam) x nypix)
    # nypix gets ceil to the next odd integer
    nypix=int(round_up_to_odd(nypix))
    ne2d=np.zeros((len(nesky),nypix))
    nefrac=np.zeros(nypix)
    for i in range(nypix):
        j=i-int(np.floor(nypix/2))
        ymin=(j-0.5)/np.sqrt(2)/(dfib/2.355)
        ymax=(j+0.5)/np.sqrt(2)/(dfib/2.355)
        nefrac[i]=(scipy.special.erf(ymax)-scipy.special.erf(ymin))/2. 
               
    for i in range(len(nesky)):
        ne2d[i,:]=nesky[i]*nefrac
        
    return ne2d
    
def mknenoise(neobj, nesky, rdn, dark, texp):
    # Return noise per pixel (generic, can do 1 pixel, a 1D or 2D spectrum depending on format of neobj and nesky)
    return np.sqrt(neobj+nesky+rdn**2+dark*texp)

def mkrandnoise(nenoise):
    return np.random.randn(*nenoise.shape)*nenoise

def round_up_to_odd(f):
    return np.ceil(f) // 2 * 2 + 1

def mkstnivar(neobj2d, nenoise2d, aper):
    # Returns S/N, S, and, N of 1D spectrum after collapsing using optimal extraction 
    stnout=np.zeros(neobj2d.shape[0])
    sout=np.zeros(neobj2d.shape[0])
    nout=np.zeros(neobj2d.shape[0])
    dy=np.arange(neobj2d.shape[1])-np.floor(neobj2d.shape[1]/2.)
    sely=np.abs(dy)<=aper/2.
    for i in range(neobj2d.shape[0]):
        faux=neobj2d[i,sely]/neobj2d[i,sely]*np.sum(neobj2d[i,sely])
        eaux=nenoise2d[i,sely]/neobj2d[i,sely]*np.sum(neobj2d[i,sely])
        sout[i]=np.sum(faux/eaux**2)/np.sum(1./eaux**2)
        nout[i]=np.sqrt(1.0/np.sum(1./eaux**2))
        stnout[i]=sout[i]/nout[i]

    return [stnout, sout, nout]
        
def mkstnaper(neobj2d, nenoise2d, aper):
    # Returns S/N, S, and, N of 1D spectrum after collapsing using sum over aperture 
    stnout=np.zeros(neobj2d.shape[0])
    sout=np.zeros(neobj2d.shape[0])
    nout=np.zeros(neobj2d.shape[0])
    dy=np.arange(neobj2d.shape[1])-np.floor(neobj2d.shape[1]/2.)
    sely=np.abs(dy)<=aper/2.
    for i in range(neobj2d.shape[0]):
        faux=neobj2d[i,sely]
        eaux=nenoise2d[i,sely]
        sout[i]=np.sum(faux)
        nout[i]=np.sqrt(np.sum(eaux**2))
        stnout[i]=sout[i]/nout[i]

    return [stnout, sout, nout]

# ==================================================
# ==================================================

# PROGRAM STARTS HERE

# LCOETC Database directory
dbdir='database'

# Temporary output directory
#tmpdir='/Users/gblanc/lcoetc/tmp'
tmpdir='output'

# Remove temporary files
if os.path.isfile(tmpdir+"/lvmetc_out.dat"):
    os.remove(tmpdir+"/lvmetc_out.dat")
if os.path.isfile(tmpdir+"/lvmetc_out.fits"):
    os.remove(tmpdir+"/lvmetc_out.fits")
if os.path.isfile(tmpdir+"/lvmetc_out.png"):
    os.remove(tmpdir+"/lvmetc_out.png")

c=2.99792458e18   # Speed of light in [A/s]

# Check time and start log file
startTime = datetime.now()
logf=open(tmpdir+'/lvmetc_log.txt', 'w')
logf.write(str(startTime)+' - Start\n')

#Deals with inputing data into python from the html form
#form = cgi.FieldStorage()


#read in the configuration parameters
configfile="config_file.txt"
d = {}
with open(configfile) as f:
    for line in f:
       (key, val) = line.split()
       d[key] = val

# Define parameters from HTML Form
ofile=dbdir+'/object/'+d["template"]
abmag=float(d["abmag"])                              # [mag arcsec2]
filfile=dbdir+'/filters/'+d["tempfilter"]

elamfile=dbdir+'/zeropoints/LVM_'+d["telescope"]+'_ELAM.dat'

binspec=int(d["binspec"])          # Binning in spectral direction [-]
binspat=int(d["binspat"])          # Binning in spatial direction [-]

nmoon=int(d["nmoon"])              # [-]
amass=float(d["amass"])            # [amass]
dpsf=float(d["dpsf"])              # [arcsec]
sfile=dbdir+'/sky/LVM_'+d["telescope"]+'_SKY_'+str(nmoon)+'.dat'

texp=float(d["texp"])              # [s]
nexp=int(d["nexp"])                # [-]
aper=float(d["aper"])              # [arcsec]

# Transmission file to use
klamfile=dbdir+'/sky/LVM_'+d["telescope"]+'_KLAM.dat'

# Read parameter file for TELESCOPE
paramfile=dbdir+'/instruments/LVM_'+d["telescope"]+'_PARAM.dat'

params=ascii.read(paramfile)

atel=float(params['col2'][params['col1']=='ATEL'])
lmin=float(params['col2'][params['col1']=='LMIN'])
lmax=float(params['col2'][params['col1']=='LMAX'])
lrange=[lmin, lmax]
dspax=float(params['col2'][params['col1']=='DSPAX'])
dfib=float(params['col2'][params['col1']=='DFIB'])
ddisp=float(params['col2'][params['col1']=='DDISP'])
dlsf=float(params['col2'][params['col1']=='DLSF'])
gain=float(params['col2'][params['col1']=='GAIN'])
rdn=float(params['col2'][params['col1']=='RDN'])
dark=float(params['col2'][params['col1']=='DARK'])

# 2D spectra are 5 times the fiber FWHM tall
nypix=dfib*5.                                 

dfiblam=dfib*ddisp

# Apply binning in spatial and spectral directions
dfib=dfib/binspat
ddisp=ddisp*binspec


logf.write(str(datetime.now() - startTime)+" - Done Defining Parameters\n")


## Create arrays 

# Wavelength
lam=mklam(lrange, ddisp)

logf.write(str(datetime.now() - startTime)+" - Done making lam array\n")

# Atmospheric Extinction Coefficient
klam=mkklam(lam, dlsf, dfiblam, klamfile)

logf.write(str(datetime.now() - startTime)+" - Done making klam array\n")

# Sky spectrum
slam=mkslam(lam, dlsf, dfiblam, sfile)*dspax**2

logf.write(str(datetime.now() - startTime)+" - Done making slam array\n")

# Object Continuum Spectrum (scaled by spaxel area)
if d['template'] != 'flat':
    olam=mkolam(lam, dlsf, dfiblam, ofile, filfile, abmag)*dspax**2
else:
    olam=10.0**(-0.4*(abmag+48.6))*c/lam**2*dspax**2

logf.write(str(datetime.now() - startTime)+" - Done making flux arrays\n")
    
    
# Object Added Emission Line
if d['addline'] == "1":
    linelam=float(d['linelam'])
    lineflux=float(d['lineflux'])
    linefwhm=float(d['linefwhm'])
    tmplam=mklinelam(lam, dlsf, dfiblam, linelam, lineflux, linefwhm)*dspax**2
    olam=olam+tmplam

logf.write(str(datetime.now() - startTime)+" - Done adding line to klam array\n")

    
# System Throughput
elam=mkelam(lam, elamfile)

logf.write(str(datetime.now() - startTime)+" - Done making elam array\n")


# Calculate electron counts and S/N for single exposure
neobj=mkneobj(lam, olam, ddisp, texp, elam, atel, klam, amass)
neobj2d=mkneobj2d(neobj, dfib, nypix)
nesky=mknesky(lam, slam, ddisp, texp, elam, atel)
nesky2d=mknesky2d(nesky, dfib, nypix)
nenoise2d=mknenoise(neobj2d, nesky2d, rdn, dark, texp)

logf.write(str(datetime.now() - startTime)+" - Done making electron arrays\n")

# The test below shows that 1.5*FWHM is ideal for APER and 2.5*FWHM for OPTIMAL. Using 2*FWHM
#for i in np.arange(1, 3, 0.1):
#    outivar=mkstnivar(neobj2d, nenoise2d, aper=i*dpsf/dpix)
#    outaper=mkstnaper(neobj2d, nenoise2d, aper=i*dpsf/dpix)
#    print(np.median(outivar[0]), np.median(outaper[0]))

#outivar=mkstnivar(neobj2d, nenoise2d, aper)
outaper=mkstnaper(neobj2d, nenoise2d, aper)
outapersky=mkstnaper(nesky2d, nenoise2d, aper)
#npixaper=int(np.median(outapersky[1]/nesky))

# Calculate electron counts and S/N for coadd of nexp exposure

neobj2dcoadd=nexp*neobj2d
nenoise2dcoadd=np.sqrt(nexp)*nenoise2d
#outivarcoadd=mkstnivar(neobj2dcoadd, nenoise2dcoadd, aper=2*dpsf/dpix)
outapercoadd=mkstnaper(neobj2dcoadd, nenoise2dcoadd, aper)

logf.write(str(datetime.now() - startTime)+" - Done making S/N arrays\n")

# Write ASCII files with output 1D spectra / noise / S/N / realization

#ascii.write([lam, outivar[0], outivar[1], outivar[2], outaper[0], outaper[1], outaper[2], outapersky[1]], tmpdir+'/lvmetc_out.dat', format='commented_header', names=['Wavelength_[A]', 'S/N_Optimal', 'S_Optimal_[e]', 'N_Optimal_[e]', 'S/N_Aperture', 'S_Aperture_[e]', 'N_Aperture_[e]', 'SKY_Aperture_[e]'], formats={'Wavelength_[A]':'%.3f', 'S/N_Optimal':'%.2f', 'S_Optimal_[e]':'%.2f', 'N_Optimal_[e]':'%.2f', 'S/N_Aperture':'%.2f', 'S_Aperture_[e]':'%.2f', 'N_Aperture_[e]':'%.2f', 'SKY_Aperture_[e]':'%.2f'})

ascii.write([lam, outaper[0], outaper[1], outaper[2], outapersky[1], outapercoadd[0], outapercoadd[1], outapercoadd[2], nexp*outapersky[1]], tmpdir+'/lvmetc_out.dat', format='commented_header', names=['Wavelength_[A]', 'S/N_Aperture', 'S_Aperture_[e]', 'N_Aperture_[e]', 'SKY_Aperture_[e]', 'S/N_Aperture_Coadd', 'S_Aperture_Coadd_[e]', 'N_Aperture_Coadd_[e]', 'SKY_Aperture_Coadd_[e]'], formats={'Wavelength_[A]':'%.3f', 'S/N_Aperture':'%.2f', 'S_Aperture_[e]':'%.2f', 'N_Aperture_[e]':'%.2f', 'SKY_Aperture_[e]':'%.2f', 'S/N_Aperture_Coadd':'%.2f', 'S_Aperture_Coadd_[e]':'%.2f', 'N_Aperture_Coadd_[e]':'%.2f', 'SKY_Aperture_Coadd_[e]':'%.2f'})

logf.write(str(datetime.now() - startTime)+" - Done writing ASCII\n")

# Write FITS files with simulated data (in units of electrons)

hdu0=fits.PrimaryHDU()
hdu1=fits.ImageHDU()
hdu2=fits.ImageHDU()
hdu3=fits.ImageHDU()
hdu4=fits.ImageHDU()
hdul=fits.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4])


hdul.writeto(tmpdir+'/lvmetc_out.fits', overwrite=1)

fits.update(tmpdir+'/lvmetc_out.fits', np.transpose(neobj2d+nesky2d+mkrandnoise(nenoise2d)), 1)
fits.update(tmpdir+'/lvmetc_out.fits', np.transpose(neobj2d), 2)
fits.update(tmpdir+'/lvmetc_out.fits', np.transpose(nesky2d), 3)
fits.update(tmpdir+'/lvmetc_out.fits', np.transpose(nenoise2d), 4)

hdul=fits.open(tmpdir+'/lvmetc_out.fits')

for i in np.arange(4):
    hdul[i+1].header.set('DISPAXIS', '1')
    hdul[i+1].header.set('CTYPE1', 'LINEAR')
    hdul[i+1].header.set('CUNIT1', 'Angstroms')
    hdul[i+1].header.set('CRPIX1', '1')
    hdul[i+1].header.set('CRVAL1', str(lam[0]))
    hdul[i+1].header.set('CD1_1', str(ddisp))


hdul.writeto(tmpdir+'/lvmetc_out.fits', overwrite=1)
hdul.close()

logf.write(str(datetime.now() - startTime)+" - Done writing FITS\n")


# Make Plots
plt.close()
plt.figure(1, figsize=(12,10))

plt.subplot(221)
#plt.plot(lam, outivar[0], 'b', label='S/N Optimal')
plt.plot(lam, outaper[0], 'r', label='S/N Aperture')
plt.legend(loc=0, fontsize=10)
plt.title('Single '+str(texp)+' sec. Exposure - '+str(aper)+' Pixel Extraction Aperture')
plt.ylabel("S/N per binned spectral pixel")
plt.xlabel('Wavelength')
#plt.axis([5000, 6000, 0, 10])

plt.subplot(223)
# plt.plot(lam, outivar[1], 'b--', label='SIGNAL Optimal')
# plt.plot(lam, outivar[2], 'g--', label='NOISE Optimal')
plt.plot(lam, outaper[1], 'b', label='SIGNAL Aperture')
plt.plot(lam, outaper[2], 'g', label='NOISE Aperture')
plt.plot(lam, outapersky[1], 'r', label='SKY Aperture')
plt.yscale('log')
plt.legend(loc=0, fontsize=10)
plt.ylabel("electrons per binned spectral pixel")
plt.xlabel('Wavelength')

plt.subplot(222)
#plt.plot(lam, outivarcoadd[0], 'b', label='S/N Optimal')
plt.plot(lam, outapercoadd[0], 'r', label='S/N Aperture')
plt.yscale('linear')
plt.legend(loc=0, fontsize=10)
plt.title(str(nexp)+' x '+str(texp)+' sec. Exposures - '+str(aper)+' Pixel Extraction Aperture')
plt.ylabel("S/N per binned spectral pixel")
plt.xlabel('Wavelength')
#plt.axis([5000, 6000, 0, 10])

plt.subplot(224)
# plt.plot(lam, outivarcoadd[1], 'b--', label='SIGNAL Optimal')
# plt.plot(lam, outivarcoadd[2], 'g--', label='NOISE Optimal')
plt.plot(lam, outapercoadd[1], 'b', label='SIGNAL Aperture')
plt.plot(lam, outapercoadd[2], 'g', label='NOISE Aperture')
plt.plot(lam, nexp*outapersky[1], 'r', label='SKY Aperture')
plt.yscale('log')
plt.legend(loc=0, fontsize=10)
plt.ylabel("electrons per binned spectral pixel")
plt.xlabel('Wavelength')

logf.write(str(datetime.now() - startTime)+" - Done Plotting\n")


plt.tight_layout()
plt.savefig(tmpdir+'/lvmetc_out.png')
plt.show()
#plt.savefig(sys.stdout, format='png')
logf.write(str(datetime.now() - startTime)+" - Done Outputing Plots\n")

#zp2=elam2zp(0.16, lam, 1.0, 1.0, 1.0, atel, klam, 1.0)


logf.write(str(datetime.now() - startTime)+" - End\n")
logf.close()

