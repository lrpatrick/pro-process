"""

PROMETEO functions to read ASCII files and create the standardised PROMETEO
FITS and ASCII files

This script should be used to generate the input for pfits 

Changelog

27-05-2020:

This file has been updated to include

1.  Functionality for the y-error to be included in the input spectrum.
getdata() has been extended to do this.  Antonio, do you know a better way to
do this?

2. The new header keywords. The list cat_kw has been changed

Author: LRP
Date: 31-03-2020
"""

import astropy.io.ascii as ascii
import astropy.io.fits as fits
import glob
import numpy as np
import os

import pfits


def getdata(data_path):
    """
    A function to get data from an input file

    There must be a better way to test if a file has two or three columns ...
    """
    try:
        # Assume file has three columns
        spec1 = ascii.read(data_path, names=['x', 'y', 'yerr'])
        data = np.column_stack((spec1['x'], spec1['y'], spec1['yerr']))
    except:
        # If not, it must have two columns
        spec1 = ascii.read(data_path, names=['x', 'y'])
        yerr = np.zeros_like(spec1['x'])
        data = np.column_stack((spec1['x'], spec1['y'], yerr))

    return data



# Get path of input catalogue
# pcat = input('[INFO] Please enter full path of input catalogue:\n')

# pcat = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/30_lmc_ccd2/30_lmc_ccd2_catalog_pro.dat'
# pcat = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/pro_input/AAO-cat.dat'
pcat = 'AAO-cat-test.csv'
# Path for spectra
# For the interface this will be taken from uploaded archive
tree = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/pro_input/'
# Required keywords to be referenced
cat_kw = ['SPECTRUM', 'OBJECT', 'DATE-OBS', 'RA', 'DEC',
          'TEL', 'INS', 'TEXP', 'ID-SIM', 'RV-FLAG',
          'RV', 'RV-ERR', 'RV-REF']
# Read input catalogue
cat = ascii.read(pcat, names=cat_kw)

for i, fascii in enumerate(cat['SPECTRUM']):
    print('[INFO] Processing file {}'.format(fascii))
    # pid_i = pid_start + i
    print('[INFO] Extracting data from input FITS')
    try:
        # Get Data:
        data = getdata(tree + cat['SPECTRUM'][i])
        # Get header information from catalogue
        keyword_values = [cat[i][kw] for kw in cat_kw[1:]]

    except:
        print('[WARNING] Unable to generate appropriate data')
    try:
        new_hdu = pfits.write_fits(data, keyword_values, oldfits=None)
    except:
        print('[WARNING] Unable to write new FITS file')
    try:
        ascii_data = pfits.write_ascii(data, keyword_values)
    except:
        print('[WARNING] Unable to write new ASCII file')


# print('[INFO] Final PID {}'.format(pid_i))

# Update database README file
# final_pid = pid_i
# readme_fin = pfits.write_proid(readme_path, final_pid)