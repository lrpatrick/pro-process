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


def getdata(x_option, y_option, yerr_option, fits):
    """
    A function to get and format data depending on the type of input

    More flexibility should be added.


    """
    # Flux
    if y_option == 0:
        y = fits[y_ext].data.field(yname)[0]
    elif y_option == 1:
        y = fits[y_ext].data[0] 

    # Wavelength
    if x_option == 0:
        x = fits[x_ext].data.field(xname)[0]
    elif x_option == 1:
        x = pfits.lam_axis(fits[crval1_ext].header[crval1],
                   fits[cdelt1_ext].header[cdelt1],
                   y.shape[0])
    # Flux error
    if yerr_option == 0:
        yerr = np.zeros_like(x)
    elif yerr_option == 1:
        yerr = np.zeros_like(x)
        # yerr = fits[yerr_ext].data.field(yerrname)[0]
    elif yerr_option == 2:
        yerr = np.zeros_like(x)
        # yerr = fits[y_ext].data[0] 

    # Stack data in useful format:
    data = np.column_stack((x, y, yerr))
    return data

# HARPS
flist = glob.glob('/home/lee/Work/ngc330/NGC330_HARPS/*.fits')
# IACOB
# flist = glob.glob('/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/IACOB/raw/*.fits')
# Cygnus
# flist = glob.glob('/home/lee/Work/PROMETEO/cygOB2_spectra/spectra_CYGOB2_Sara/INT/spn_ag16/*.fits')
pid_start = int(input('[INFO] Please enter ID number for first spectrum as integer:\n')) or int(0)

# User input KEYWORDS for all spectra:

# TODO: Guess ...

# The HARPS keywords will be the default values
# Input keywords 
# HARPS
in_kws = ['OBJECT', 'DATE-OBS', 'RA', 'DEC',
          'TELESCOP', 'INSTRUME', 'EXPTIME']
# Input extensions
in_exts = np.zeros(len(in_kws), dtype=int)

# IACOB
# in_kws = ['OBJECT', 'DATE-OBS', 'I-RA', 'I-DEC',
#           'TELESCOP', 'INSTRUME', 'I-TEXP']

# Cygnus
# in_kws = ['OBJECT', 'DATE-OBS', 'RA', 'DEC',
#           'TELESCOP', 'INSTRUME', 'TEXP']

# Input extensions
in_exts = np.zeros(len(in_kws), dtype=int)

# Input data
print('[INFO] Data structure definition:')
print('[INFO] Wavelength:')
x_option = int(input('[INPUT] Options:\n 0 : Array Name\n 1 : Keywords\n')) or 0
if x_option == 0:
    xname = input('Array name:\n')
    x_ext = int(input('Header extension:\n'))
elif x_option == 1:
    crval1 = input('Reference Position header keyword:\n')
    crval1_ext = int(input('Reference Position header keyword extension:\n'))
    cdelt1 = input('Pixel scale keyword:\n')
    cdelt1_ext = int(input('Pixel scale keyword extension:\n'))
else:
    print('[WARNING] INPUT not understood')

print('[INFO] Flux:')
y_option = int(input('[INPUT] Options:\n 0 : Array Name\n 1 : Extension\n')) or 0
if y_option == 0:
    yname = input('Array name:\n')
    y_ext = int(input('Header extension:\n'))
elif y_option == 1:
    y_ext = int(input('Header extension:\n'))
else:
    print('[WARNING] INPUT not understood')

print('[INFO] Flux Error:')
yerr_option = int(input('[INPUT] Options:\n 0 : Zeros\n 1 : Array Name\n 2: Extension\n')) or 0
if yerr_option == 0:
    pass
elif yerr_option == 1:
    yname = input('Array name:\n')
    y_ext = int(input('Header extension:\n'))
elif yerr_option == 2:
    y_ext = int(input('Header extension:\n'))
else:
    print('[WARNING] INPUT not understood')


for i, ffits in enumerate(flist[:20]):
    print('[INFO] Processing file {}'.format(ffits))
    pid_i = pid_start + i
    print('[INFO] Extracting data from input FITS')
    try:
        infits = fits.open(ffits)
        # Get Data:
        data = getdata(x_option, y_option, yerr_option, infits)
        # Get header information from in_kws
        keyword_values = [infits[ext].header[kw]
                          for ext, kw in zip(in_exts, in_kws)]

    except:
        print('[WARNING] Unable to generate appropriate data')

    try:
        new_hdu = pfits.write_fits(pid_i, data, keyword_values, infits)
    except:
        print('[WARNING] Unable to write new FITS file')
    try:
        ascii_data = pfits.write_ascii(pid_i, data, keyword_values)
    except:
        print('[WARNING] Unable to write new ASCII file')

    # print('[INFO] Write output ascii file:')


print('[INFO] Final PID {}'.format(pid_i))
