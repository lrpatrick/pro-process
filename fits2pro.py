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


def getdata(data_options, data_keywords, fits):
    """
    A function to get and format data depending on the type of input

    More flexibility should be added.


    """
    x_option, y_option, yerr_option = data_options

    # Flux
    if y_option == 0:
        y = fits[data_keywords[5]].data.field(data_keywords[4])[0]
    elif y_option == 1:
        y = fits[data_keywords[5]].data[0] 

    # Wavelength
    if x_option == 0:
        x = fits[data_keywords[1]].data.field(data_keywords[0])[0]
    elif x_option == 1:
        x = pfits.lam_axis(fits[data_keywords[1]].header[data_keywords[0]],
                   fits[data_keywords[3]].header[data_keywords[2]],
                   y.shape[0])
    # Flux error
    if yerr_option == 0:
        yerr = np.zeros_like(x)
    elif yerr_option == 1:
        yerr = np.zeros_like(x)
        # yerr = fits[yerr_ext].data.field(yerrname)[0]
    elif yerr_option == 2:
        yerr = np.zeros_like(x)
        # yerr = fits[data_keywords[5]].data[0] 

    # Stack data in useful format:
    data = np.column_stack((x, y, yerr))
    return data


def guess_fits(infits):
    """
    Function to guess the structure of a FITS file
    
    What all the fits have in common is: 
        0th extension 
        infits[0].header -- in every fits
        infits[0].dat    -- in every fits, but may be None
    infits[0].header
        All FITS have an EXTEND keyword (or they should) -- but it is
        sometimes unreliable!
        Every fits has some sort of NAXIS keyword that gives us an indication
        of what type of FITS it is
        NAXIS == 0 indicates that there is no data in the 0th extension
        (this is part of a standard followed by the Virtual Observatory)
    
    Arguments: 
    infits : astropy.io.fits.hdu.hdulist.HDUList
        One FITS file that represents of all of the files in the archive
    
    Returns:
    data_options : list
        A list of options to describe the input data. This consists of three
        integers corresponding to: x_option, y_option, yerr_option

    data_keywords : list
        A list of header keywords and names corresponding to the data choices.

    header_kws : list
        A list of header keywords that will be used to get header information.
        This list contains six strings with the default keywords.

    header_exts : list
        A list of header extensions. The default is that these are all zero.

    """

    data_options, data_keywords, header_kws, header_exts = pfits.interface_defaults()
    data_options, data_keywords = pfits.get_data_structure(infits)
    header_kws, header_exts = pfits.get_header_structure(infits)
    return data_options, data_keywords, header_kws, header_exts


# HARPS
flist = glob.glob('/home/lee/Work/ngc330/NGC330_HARPS/*.fits')
# IACOB
# flist = glob.glob('/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/IACOB/raw/*.fits')
# Cygnus
# flist = glob.glob('/home/lee/Work/PROMETEO/cygOB2_spectra/spectra_CYGOB2_Sara/INT/spn_ag16/*.fits')

# Get ID number:
# pid_start = int(input('[INFO] Please enter ID number for first spectrum as integer:\n')) or int(0)
readme_path = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/ASTRO+database/README'
pid_start = pfits.get_proid(readme_path)

# User input KEYWORDS for all spectra:

# Guess input data structure

data_options, data_keywords,\
    header_kws, header_exts = guess_fits(fits.open(flist[0]))

# The HARPS keywords will be the default values
# Input keywords 
# HARPS
# in_kws = ['OBJECT', 'DATE-OBS', 'RA', 'DEC',
#           'TELESCOP', 'INSTRUME', 'EXPTIME']
# # Input extensions
# in_exts = np.zeros(len(in_kws), dtype=int)

# # IACOB
# # in_kws = ['OBJECT', 'DATE-OBS', 'I-RA', 'I-DEC',
# #           'TELESCOP', 'INSTRUME', 'I-TEXP']

# # Cygnus
# # in_kws = ['OBJECT', 'DATE-OBS', 'RA', 'DEC',
# #           'TELESCOP', 'INSTRUME', 'TEXP']

# # Input extensions
# in_exts = np.zeros(len(in_kws), dtype=int)

# Input data
print('[INFO] Data structure definition:')
print('[INFO] Wavelength:')
print('[INFO] From guessing function x-option is {0} with the following option\n {1} ext:{2} {3} ext:{4}'
    .format(data_options[0], data_keywords[0], data_keywords[1], data_keywords[2], data_keywords[3]))
x_option = input('[INPUT] Update the above choice? Options:\n 0 : Array Name\n 1 : Keywords\nHit return to use above option') or None
if x_option == None:
    x_option = data_options[0]
else:
    x_option = int(x_option)
    data_options[0] = x_option
    if x_option == 0:
        data_keywords[0] = input('Array name:\n')
        data_keywords[1] = int(input('Header extension:\n'))
        data_keywords[2] = None
        data_keywords[3] = None
    elif x_option == 1:
        data_keywords[0] = input('Reference Position header keyword:\n')
        data_keywords[1] = int(input('Reference Position header keyword extension:\n'))
        data_keywords[2] = input('Pixel scale keyword:\n')
        data_keywords[3] = int(input('Pixel scale keyword extension:\n'))
    else:
        print('[WARNING] INPUT not understood')

print('[INFO] Flux:')
print('[INFO] From guessing function y-option is {0} with the following option\n {1} ext:{2}'
    .format(data_options[1], data_keywords[4], data_keywords[5]))
y_option = input('[INPUT] Update the above choice? Options:\n 0 : Array Name\n 1 : Extension\nHit return to use above option') or None
if y_option == None:
    y_option = data_options[1]
else:
    y_option = int(y_option)
    data_options[1] = y_option
    if y_option == 0:
        data_keywords[4] = input('Array name:\n')
        data_keywords[5] = int(input('Header extension:\n'))

    elif y_option == 1:
        data_keywords[5] = int(input('Header extension:\n'))        
    else:
        print('[WARNING] INPUT not understood')

print('[INFO] Flux Error:')
print('[INFO] From guessing function yerr-option is {0}'
    .format(data_options[2]))
yerr_option = input('[INPUT] Update the above choice? Options:\n 0 : Zeros\n 1 : Array Name\n 2: Extension\n') or None
if yerr_option == None:
    yerr_option = data_options[2]
else:
    if yerr_option == 0:
        pass
    elif yerr_option == 1:
        yname = input('Array name:\n')
        y_ext = int(input('Header extension:\n'))
    elif yerr_option == 2:
        y_ext = int(input('Header extension:\n'))
    else:
        print('[WARNING] INPUT not understood')

# Check Header structure
print('[INFO] Current header values:')

for i, kw in enumerate(header_kws):
    print('KW: "{}" EXT {}'.format(kw, header_exts[i]))
    update = input('[INPUT] Update the above choice? [y/N]\n') or 'n'
    if update.lower() == 'y':
        header_kws[i] = input('[INPUT] New Keyword (str):\n')
        header_exts[i] = input('[INPUT] New Extension (int):\n')


# a = input('[INFO] This is a chance to stop before entering the unbreakable loop!')


for i, ffits in enumerate(flist):
    print('[INFO] Processing file {}'.format(ffits))
    pid_i = pid_start + i
    print('[INFO] Extracting data from input FITS')
    try:
        infits = fits.open(ffits)
        # Get Data:
        data = getdata(data_options, data_keywords, infits)
        # Get header information from in_kws
        keyword_values = [infits[ext].header[kw]
                          for ext, kw in zip(header_exts, header_kws)]

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

# Update database README file
final_pid = pid_i
readme_fin = pfits.write_proid(readme_path, final_pid)


