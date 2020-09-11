import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy
from astropy.io import ascii
from astropy.table import Table
from scipy.interpolate import interp1d


for i in range(15):
    data0=ascii.read('skyflam_'+str(i)+'.dat')

    lam0=np.array(data0['col1'])
    slam0=np.array(data0['col2'])

    

# MAGE
    lmin=3300
    lmax=10100
    ddisp=0.35
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_MAGE_ECHELLETTE_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})


# LDSS3 VPHALL
    lmin=3650
    lmax=10100
    ddisp=1.89
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN2_LDSS3_VPHALL_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

# LDSS3 VPHBLUE
    lmin=3850
    lmax=6150
    ddisp=0.682
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN2_LDSS3_VPHBLUE_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

# LDSS3 VPHRED
    lmin=6000
    lmax=10100
    ddisp=1.175
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN2_LDSS3_VPHRED_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})
    
# MIKE BLUE
    lmin=3350
    lmax=5000
    ddisp=0.02
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN2_MIKE_BLUE_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})
    
    # MIKE RED
    lmin=4900
    lmax=10000
    ddisp=0.05
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN2_MIKE_RED_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})
    
    
    # IMACS F2_150_11
    lmin=4600
    lmax=9000
    ddisp=2.630
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F2_150_11_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

    
# IMACS F2_200_15
    lmin=4300
    lmax=9400
    ddisp=2.037
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F2_200_15_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})
    

    # IMACS F2_300_17
    lmin=4000
    lmax=9000
    ddisp=1.341
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F2_300_17_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})


    # IMACS F2_300_26
    lmin=4700
    lmax=9400
    ddisp=1.25
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F2_300_26_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})


    # IMACS F4_150-3_3.4
    lmin=4000
    lmax=9400
    ddisp=1.453
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F4_150-3_3.4_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

    # IMACS F4_300-4_6.0
    lmin=3600
    lmax=8000
    ddisp=0.743
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F4_300-4_6.0_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})


    # IMACS F4_600-8_9.3
    lmin=3600
    lmax=6700
    ddisp=0.378
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F4_600-8_9.3_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})


    # IMACS F4_600-13_14.0
    lmin=6700
    lmax=9000
    ddisp=0.387
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F4_600-13_14.0_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

    # IMACS F4_1200-17_19.0
    lmin=4200
    lmax=5900
    ddisp=0.194
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F4_1200-17_19.0_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

    # IMACS F4_1200-27_27.0
    lmin=6300
    lmax=7800
    ddisp=0.096
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F4_1200-27_27.0_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

    # IMACS F4_1200-27_33.5
    lmin=7900
    lmax=9300
    ddisp=0.188
    lam1=np.arange(lmin-100, lmax+ddisp+100, ddisp/3)
    f=interp1d(lam0, slam0, fill_value='extrapolate')
    slam1=f(lam1)
    data1=Table([lam1, slam1], names=('col1', 'col2'))
    ascii.write(data1, 'MAGELLAN1_IMACS_F4_1200-27_33.5_SKY_'+str(i)+'.dat', format='no_header', formats={'col1':'%.3f', 'col2':'%.3e'})

    
