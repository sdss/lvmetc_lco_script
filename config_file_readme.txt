Full set of options for configuration file inputs.  Default values are indicated below and included as a ready to use configuration file (config_file_default.txt).

Continuum Spectral Template
template
O5V_Pickles.dat         O5V star
B0V_Pickles.dat         B0V star
A0V_Pickles.dat         A0V star
F0V_Pickles.dat         F0V star (default value)
G0V_Pickles.dat         G0V star
K0V_Pickles.dat         K0V star
M0V_Pickles.dat         M0V star
flat                    Flat Continuum (AB)


Continuum Surface Brightness [AB mag/arcsec^2]]
abmag
23.5        (default value)

Filter
tempfilter
sdss_u.dat          SDSS u'
sdss_g.dat          SDSS g'
sdss_r.dat          SDSS r' (default value)
sdss_i.dat          SDSS i'
sdss_z.dat          SDSS z'

Add Emission Line
addline
0       No
1       Yes (default value)


Line Wavelenght [A]
linelam
6563.0  (default value)

Line Surface Brightness [erg/(s cm^2 arcsec^2)]
lineflux
3e-18   (default value)

Line FWHM [A]
linefwhm
5.0     (default value)

===========
Instrumental Setup
===========

Telescope
telescope
APO25       APO 2.5m
LCO25       LCO 2.5m
NMSU10      NMSU 1m
LVM160      LVM-160 0.16m (default value)

Binning Spatial Direction
binspat
1       (default value)
2
3
4
5
6
7
8

Binning Spectral Direction
binspec
1       (default value)
2
3
4
5
6
7
8

============
Observing Conditions
===========


Days from New Moon
nmoon
0       (default value)
1
2
3
4
5
6
7
8
9
10
11
12
13
14

Airmass
amass
1.2     (default value)

Seeing FWHM [arcsec] (not used)
dpsf
0.6     (default value)

==========
Exposure and Extraction Parameters
==========

Single Frame Exposure Time [seconds]
texp
300     (default value)

Number of Exposures to Coadd
nexp
3       (default value)

Fiber Extraction Aperture [N pixels]
aper
4       (default value)
