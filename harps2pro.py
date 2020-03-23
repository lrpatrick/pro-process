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

rob74_ids = {1.: 'NGC330-R74-A6', 2.: 'NGC330-R74-A3', 3.: 'NGC330-R74-A14',
             4.: 'NGC330-R74-A9', 5.: 'NGC330-R74-B31', 6.: 'NGC330-R74-B40',
             7.: 'NGC330-R74-B15', 8.: 'NGC330-R74-B42', 9.: 'NGC330-R74-A42'}


flist_pro = glob.glob('/home/lee/Work/ngc330/NGC330_HARPS/*.fits')

pid_start = int(input('[INFO] Please enter ID number for first spectrum as integer:\n')) or int(0)

for i, fharps in enumerate(flist_pro[:10]):
    print('[INFO] Processing file {}'.format(fharps))
    pid_i = pid_start + i
    print('[INFO] Write output ascii file:')
    try:
        pro1 = fits.open(fharps)
        ngc330_id = rob74_ids[int(pro1[0].header['OBJECT'].split('-')[1])]
        pro1[0].header['OBJECT'] = ngc330_id
        h1, obj_name_str = pfits.write_ascii_header_harps(pro1, pid_i)
        pro_outname_ascii = obj_name_str + '_' + os.path.basename(fharps)[:-5] + '_pro.dat'
        x = pro1[1].data.field('WAVE')[0]
        y = pro1[1].data.field('FLUX')[0]
        ascii_data = pfits.write_ascii(pro_outname_ascii, x, y,
                                       yerr=0, head=h1)
        
        print('[INFO] Write output FITS file:')
        # Start with a fresh fits header based on the data
        # Add important keywords:
        new_hdu = pfits.write_fits_harps(pro1, ascii_data, pid_i)

    except:
        print('[WARNING] FITS file likely not in standard format')



print('[INFO] Final PID {}'.format(pid_i))
