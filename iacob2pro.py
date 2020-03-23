"""

PROMETEO functions to read input fits/ASCII files and create the standardised PROMETEO
FITS files


This script creates new normalised FITS headers using the pfits routines


Author: LRP
Date: 13-02-2020
"""

import astropy.io.fits as fits
import glob
import numpy as np
import os

import pfits



flist_iacob = glob.glob('/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/IACOB/raw/*.fits')


pid_start = int(input('[INFO] Please enter ID number for first spectrum as integer:\n')) or int(0)

for i, fiacob in enumerate(flist_iacob[:20]):
    print('[INFO] Processing file {}'.format(fiacob))
    pid_i = pid_start + i
    print('[INFO] Write output ascii file:')
    try:
        iacob1 = fits.open(fiacob)
        h1 = pfits.write_ascii_header(iacob1, pid_i)
        iacob_outname_ascii = os.path.basename(fiacob)[:-5] + '_pro.dat'
        x = pfits.lam_axis(iacob1[0].header['CRVAL1'],
                           iacob1[0].header['CDELT1'],
                           iacob1[0].data[0].shape[0])
        ascii_data = pfits.write_ascii(iacob_outname_ascii, x, iacob1[0].data[0],
                                       yerr=iacob1[0].data[1], head=h1)
        
        print('[INFO] Write output FITS file:')
        # iacob_outname_fits = os.path.basename(fiacob)[:-5] + '_pro.fits'
        iacob_outname_fits = 'this_should_not_be_written.NOT'
        # Start with a fresh fits header based on the data
        # Add important keywords:
        new_hdu = pfits.write_fits(iacob1, ascii_data, pid_i)

    except:
        print('[WARNING] FITS file likely not in standard format')



print('[INFO] Final PID {}'.format(pid_i))
