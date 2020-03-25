"""

PROMETEO functions to read FITS files and create the standardised PROMETEO
FITS and ASCII files

This script should be used to generate the input for pfits 

TODO: Adapt this script to use as input the interface


Author: LRP
Date: 25-03-2020
"""

import astropy.io.fits as fits
import glob
import numpy as np
import os

import pfits

def getdata(option, fits):
    """
    A function to get and format data depending on the type of input

    More options should be added as they are encountered
    """
    if option == 0:
        x = fits[1].data.field('WAVE')[0]
        y = fits[1].data.field('FLUX')[0]
        yerr = np.zeros_like(x)
    elif option == 1:
        x = pfits.lam_axis(fits[0].header['CRVAL1'],
                   fits[0].header['CDELT1'],
                   fits[0].data[0].shape[0])
        y = fits[0].data[0]
        yerr = fits[0].data[1]

    data = np.column_stack((x, y, yerr))
    return data

# HARPS
# flist = glob.glob('/home/lee/Work/ngc330/NGC330_HARPS/*.fits')
# IACOB
flist = glob.glob('/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/IACOB/raw/*.fits')

pid_start = int(input('[INFO] Please enter ID number for first spectrum as integer:\n')) or int(0)

# User input KEYWORDS for all spectra:

# This information will be first guessed, then the user will have the final
# say on what actually get processed 

# These willbe the default values
# Input keywords 
# HARPS
in_kws = ['OBJECT', 'DATE-OBS', 'RA', 'DEC',
          'TELESCOP', 'INSTRUME', 'EXPTIME']
# Input extensions
in_exts = np.zeros(len(in_kws), dtype=int)

# IACOB
in_kws = ['OBJECT', 'DATE-OBS', 'I-RA', 'I-DEC',
          'TELESCOP', 'INSTRUME', 'I-TEXP']

# Input extensions
in_exts = np.zeros(len(in_kws), dtype=int)

# Input data
print('[INFO] How is the input data stored?')
print('[INFO] Options:\n 0 : wavelength and spectrum as separate arrays in 1st extension\n')
print('[INFO] 1 : Spectrum in 0th extension, wavelength from header keywords (CRVAL1, CDELT1)\n')
data_store = int(input('[INPUT] Please select an option')) or 0

for i, ffits in enumerate(flist[:20]):
    print('[INFO] Processing file {}'.format(ffits))
    pid_i = pid_start + i
    print('[INFO] Extracting data from input FITS')
    try:
        infits = fits.open(ffits)
        # Get Data:
        data = getdata(data_store, infits)
        # Get header information from in_kws
        keyword_values = [infits[ext].header[kw]
                          for ext, kw in zip(in_exts, in_kws)]

    except:
        print('[WARNING] Unable to generate appropriate data')

    try:
        new_hdu = pfits.write_fits(infits, pid_i, data, keyword_values)
    except:
        print('[WARNING] Unable to write new FITS file')
    try:
        ascii_data = pfits.write_ascii(pid_i, data, keyword_values)
    except:
        print('[WARNING] Unable to write new ASCII file')

    # print('[INFO] Write output ascii file:')


print('[INFO] Final PID {}'.format(pid_i))
