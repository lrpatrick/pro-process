"""

PROMETEO functions to read input data files and create the standardised
PROMETEO FITS and ASCII files


The project constsis of two main functions:

write_fits()
write_ascii()

These two functions are called from an input script and are used to write the
output files

Changelog:

2020-05-27 : 
    Added functionality for new keywords including optional RV
        Functions changed:
            interface_defaults()
                    The list "header_kws" has been changed
            write_fits_header() 
                    Lines 407-412 have been added
            write_ascii_header()
                    Now includes an if statement starting on line 540
        Code changed: ascii2pro.py and fits2pro.py

    Yerr functionality
        Functions changed:
            interface_defaults()
                The list "data_keywords" list has been changed
        Code changed: ascii2pro.py and fits2pro.py

    Output FITS structure changed
        Function changed:
            write_fits()
                Lines 442-446 have been added
    Added function get_date() to make sure date is in correct format
        Functions changed: write_fits_header() write_ascii_header()

Author: LRP
Date: 12-02-2020
"""

import datetime
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u 


def lam_axis(lmin, step, size):
    """
        Create an wavelength axis based on minimum wavelength,
        resolution element and number of elements in an array

        Arguments:
        lmin : float
            Minimum wavelength
        step : float
            Resolution element
        size : int
            size of array
    """
    lam = np.linspace(lmin, lmin + (step * (size - 1)), size)
    return lam


def interface_defaults():
    """
    Return deafult options for all input values

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
    # Data: 
    # Wavelength
    x_option = 0
    xname ='WAVE'
    x_ext = 0
    # Flux
    y_option = 0
    yname ='FLUX'
    y_ext = 0
    # Flux Error
    yerr_option = 0

    # Put together
    data_options = [x_option, y_option, yerr_option]
    data_keywords = [xname, x_ext, None, None,
                     yname, y_ext,
                     None, None]

    # Header:
    header_kws = ['OBJECT', 'DATE-OBS', 'RA', 'DEC',
                  'TELESCOP', 'INSTRUME', 'EXPTIME',
                  'ID-SIM', 'RV-FLAG', 'RV', 'RV-ERR', 'RV-REF']
    header_exts = np.zeros(len(header_kws), dtype=int)
    return data_options, data_keywords, header_kws, header_exts



def get_data_structure(infits):
    """
    Function to guess the structure of the input data from a given fits file

    This function returns all of the data options requested along with their
    associated keywords. How will this work when they have different keywords?

    Arguments: 
    infits : astropy.io.fits.hdu.hdulist.HDUList
        One FITS file that represents of all of the files in the archive
    

    Returns:
    data_options : list
        A list of options to describe the input data. This consists of three
        integers corresponding to: x_option, y_option, yerr_option

    data_keywords : list
        A list of header keywords and names corresponding to the data choices.
        Note that the length and content of this list changes depending on the
        choices

    """
    data_options, data_keywords, \
        header_kws, header_exts = interface_defaults()

    hdr_keys0 = list(infits[0].header.keys())
    exts = len(infits)

    # Wavelength:
    # CRVAL is always present if x_option == 1 is to be generated:
    crval_matches = [key for key in hdr_keys0 if 'CRVAL' in key]
    if len(crval_matches) == 1:
        print('[INFO] CRVAL keyword present')
        data_options[0] = 1
        data_keywords[0] = crval_matches[0]
        # Extension
        data_keywords[1] = 0
        cdelt_matches = [key for key in hdr_keys0 if key.startswith('CD')]
        if len(cdelt_matches) == 1:
            print('[INFO] CDELT keyword present')
            data_keywords[2] = cdelt_matches[0]
            # Extension
            data_keywords[3] = 0
        elif len(cdelt_matches) > 1:
            print('[INFO] Multiple CDELT keywords present')
            print('[INFO] Searching for exact matches')
            cdelt_exact = [key for key in cdelt_matches if 'CDELT1' in key]
            if len(cdelt_exact) == 1:
                data_keywords[2] = cdelt_exact[0]
                # Extension
                data_keywords[3] = 0
            else:
                # Try CD1_1
                cdelt_exact = [key for key in cdelt_matches if 'CD1_1' in key]
                data_keywords[2] = cdelt_exact[0]
                # Extension
                data_keywords[3] = 0
    elif len(crval_matches) == 0:
        print('[INFO] No CRVAL keyword present')
        data_options[0] = 0
    else:
        print('[INFO] Ambigious CRVAL keyword, returning default')

    # If no wavelength keywords in primary header and NAXIS is zero
    # There is nothing more of interest in this header
    naxis = infits[0].header['NAXIS']
    if naxis == 0:
        print('[INFO] Info from primary header exhausted. Trying others:')
        for ext in range(1, exts):
            hdr_keys = list(infits[ext].header.keys())
            # Directly search for header names
            array_names = [infits[ext].header[key] for key in hdr_keys if "TTYPE" in key]
            if array_names[0] == 'WAVE' and array_names[1] == 'FLUX':
                print('[INFO] Found a WAVE and FLUX array')
                data_options[0] = 0
                data_keywords[0] = 'WAVE'
                data_keywords[1] = ext
                data_options[1] = 0
                data_keywords[4] = 'FLUX'
                data_keywords[5] = ext
                break
            else:
                print('[INFO] No suitable keywords found. Using defaults')
    # Flux
    if infits[0].data is not None:
        print('[INFO] Data in primary header')
        data_options[1] = 1
        data_keywords[5] = 0
    else:
        print('[INFO] Data information not understood. Return defaults')
    return data_options, data_keywords



def get_header_structure(infits):
    """
    Function to retrun guessess of all the required keywords based on the
    input FITS file

    Arguments: 
    infits : astropy.io.fits.hdu.hdulist.HDUList

        One FITS file that is representative of all of the files in the
        dataset. 

    Returns:
    header_kws : list
        A list of header keywords that will be used to get header information.
        This list contains six strings with the default keywords.

    header_exts : list
        A list of header extensions. The default is that these are all zero.

    """
    # Start with the default options:
    data_options, data_keywords, \
        header_kws, header_exts = interface_defaults()

    exts = len(infits)
    for i, kw in enumerate(header_kws):
        print(kw)
        ext = 0
        # List all header keywords
        hdr_keys0 = list(infits[0].header.keys())
        # How many extensions?

        # Search for exact matches to default keywords:
        matches = [obj for obj in hdr_keys0 if kw in obj]
        # Keyword hit:

        if len(matches) == 1:
            print('[INFO] Unique keyword hit. Using this keyword and extension')
            header_kws[i] = matches[0]
            header_exts[i] = ext
        elif len(matches) > 1:
            print('[INFO] Multiple matches found')
            # Search for IACOB match
            i_match = [obj for obj in matches if 'I-' + kw == obj]
            if len(i_match) == 1:
                print('[INFO] IACOB keyword found')
                header_kws[i] = i_match[0]
                header_exts[i] = ext
            # Search for exact matches
            exact_match = [obj for obj in matches if kw == obj]
            if len(exact_match) == 1:
                print('[INFO] Exact keyword match found')
                header_kws[i] = exact_match[0]
                header_exts[i] = ext
        else:
            print('[INFO] No suitable matches')

    return header_kws, header_exts



def get_proid(readme_path):
    """
    Function to get the next ID number to be written to the files
    
    Arguments:
    readme_path : string
        Full path to the readme file to be read

    Returns:
    pid_start : int
        ID number extractedfrom input readme file
    """
    readme = ascii.read(readme_path)
    pid_start = int(readme['col2'][0])
    return pid_start


def write_proid(readme_path, pid_fin):
    """
    Funtion to write the ID number to a README file
    
    Arguments:
    readme_path : string
        Full path to the readme file to be read

    pid_fin : int
        ID number to be written to readme file readme file

    Returns:
    readme : astropy.table.table.Table
        Output table that has been written to file
    """
    readme = ascii.read(readme_path)
    readme['col2'][0] = int(pid_fin)
    ascii.write(readme, readme_path, overwrite=True)
    return readme


def write_filename(ext, name, date_obs, tel, ins, spec_start, spec_end):
    """
    Function to create the filename for the processed files

    Arguments:

    ext : str
        extension name e.g. '.fits' or '.dat'

    name : str
        Name of object. This name should not contain an underscore '_'

    date_obs : str
        Time and Date of the observation in the format 2010-09-29T00:33:22.070

    tel : str
        Telescope name

    ins : str
        Instrument name

    spec_start : str
        The start of the spectral range (in Angstroms only)
    spec_end : str
        The end of the spectral range (in Angstroms only)

    Returns: 

    fname : str
        The final standardised filename using only the input information

    """
    # Remove any spaces in name
    name = name.replace(" ", "")
    fname = str(name) + '_' + str(date_obs) + '_' + str(tel) + '_' + str(ins)\
            + '_' + str(spec_start) + '_' + str(spec_end) + ext

    return fname


def write_fits_header(newfits, keywords):
    """
    Create a fits header and instert it into the new fits fits
    
    Arguments:

    newfits : astropy.io.fits.hdu.hdulist.HDUList
        A fits file for the keywords to be added

    keywords : list
        A list of keywords values in following order:

        OBJECT   : str
            Name of object
        DATE-OBS : str
            Date of observation       [YYYY-MM-DDTHH:MM:SS.sss]
        RA       : float
            Right ascension of target [DEGREE]
        DEC      : float
            Declination of target     [DEGREE]
        TELESCOPE : str
            Telescope used            [free text]
        INSTRUMET : str
            Instrument used           [free text]
        EXPTIME   : float
            Exposure time             [Seconds]
        RV-FLAG   : int
            Flag to determine whether or not RV is given             
        RV        : float (optional)
            Input RV measurement
        RV-ERR    : float (optional)
            Input RV measurement error
        RV-REF    : str (optional)
            Input RV measurement reference

    Returns:

    newfits : astropy.io.fits.hdu.hdulist.HDUList
        The input newfits file with the header keywords added

    """
    # Generate date-time to put in header
    # date = str(datetime.date.today())
    date = get_date()
    comment = 'ASTRO+ processed FITS file'
    # Put data into header
    newfits.header['DATE'] = date
    # newfits.header['PROID'] = pid
    newfits.header['OBJECT'] = keywords[0]
    newfits.header['DATE-OBS'] = keywords[1]
    newfits.header['RA'] = keywords[2]
    newfits.header['DEC'] = keywords[3]
    newfits.header['TEL'] = keywords[4]
    newfits.header['INS'] = keywords[5]
    newfits.header['TEXP'] = keywords[6]
    newfits.header['ID-SIM'] = keywords[7]
    newfits.header['RV-FLAG'] = keywords[8]
    if keywords[8] == 1:
        newfits.header['RV'] = keywords[9]
        newfits.header['RV-ERR'] = keywords[10]
        newfits.header['RV-REF'] = keywords[11]
    newfits.header['COMMENT'] = comment

    return newfits


def write_fits(data, keywords, oldfits=None):

    """
    Function to write the standard PROMETEO FITS files

    Arguments:

    data : numpy.ndarray
        Array containing three columns:
            x : wavelength 
            y : flux 
            yerr : Uncertainity in flux

    keywords : list
        A list of keywords values
    
    oldfits : astropy.io.fits.hdu.hdulist.HDUList
        A fits file to be used to propagate header information

    Returns:
    prohdu : astropy.io.fits.hdu.hdulist.HDUList
        New FITS file
    """
    # Binary table data:
    new_hdu = fits.PrimaryHDU()
    bin_tab = fits.BinTableHDU.from_columns(
        [fits.Column(name='WAVE', format='E', array=data[:, 0]),
         fits.Column(name='FLUX', format='E', array=data[:, 1]),
         fits.Column(name='ERR', format='E', array=data[:, 2])])

    # Add header information to new fits
    new_hdu = write_fits_header(new_hdu, keywords)
    # Get the filename to be written out
    fname = write_filename('.fits', keywords[0], keywords[1],
                                    keywords[4], keywords[5],
                                    int(round(data[:, 0][0])),
                                    int(round(data[:, 0][-1])))
    # Add additional information from FITS if appropriate
    if oldfits is not None:
        image_hdu = fits.ImageHDU(None, header=oldfits[0].header)
    else:
        image_hdu = fits.ImageHDU(None)

    # Combine extensions into one file
    prohdu = fits.HDUList([new_hdu, bin_tab, image_hdu])

    # Write output FITS file
    prohdu.writeto('data/' + fname, output_verify='ignore')
    print('[INFO] FITS file written to {0}'.format(fname))
    return prohdu


def write_ascii(data, keywords):
    """
    Write the standard PROMETEO ascii files

    Arguments

    data : numpy.ndarray
        A three column data array containing the following:
        x : numpy.ndarray
            Continually increasing wavelength axis for spectroscopic data
        y : numpy.ndarray
            Normalised spectroscopic data. Same size as x
        yerr : numpy.ndarray
            Uncertainity on y. Same size as x

    keywords : list
        A list of keywords values

    Return 

    data : numpy.ndarray
        The unchanged input data array

    """
    fname = write_filename('.dat', keywords[0], keywords[1],
                                    keywords[4], keywords[5],
                                    int(round(data[:, 0][0])),
                                    int(round(data[:, 0][-1])))
    head1 = write_ascii_header(keywords)
    np.savetxt('data/' + fname, data, header=head1, fmt='%8.8f')
    print('[INFO] ASCII file written to {0}'.format(fname))
    return data


def write_ascii_header(keywords):
    """
    Write ASCII standard header format
    
    Arguments:

    keywords : list
        A list of keywords values in following order:

        OBJECT : str
            Name of object
        DATE-OBS  : str
            Date of observation       [YYYY-MM-DDTHH:MM:SS.sss]
        RA : float
            Right ascension of target [DEGREE]
        DEC : float
            Declination of target     [DEGREE]
        TELESCOPE : str
            Telescope used            [free text]
        INSTRUMET : str
            Instrument used           [free text]
        EXPTIME   : float
            Exposure time             [Seconds]
        RV-FLAG   : int
            Flag to determine whether or not RV is given   

    Returns:

    hascii : str
        ASCII header to be written to file

    """ 
    # Generate data to write
    date = get_date()
    # pid = str(pid)
    # Determine if RV provided
    if keywords[8] == 0:
        hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + date \
            + '\nOBJECT: ' + str(keywords[0]) \
            + '\nDATE-OBS: ' + str(keywords[1]) \
            + '\nRA: ' + str(keywords[2]) \
            + '\nDEC: ' + str(keywords[3]) \
            + '\nTELESCOPE: ' + str(keywords[4]) \
            + '\nINSTRUMENT: ' + str(keywords[5]) \
            + '\nEXPOSURE TIME: ' + str(keywords[6]) \
            + '\nID-SIM: ' + str(keywords[7]) \
            + '\nRV-FLAG: ' + str(keywords[8]) \
            + '\nCOL1: WAVELENGH ARRAY [Angstrom]' \
            + '\nCOL2: NORMALISED FLUX ARRAY [arbitary units]' \
            + '\nCOL3: ERROR SPECTRUM [arbitary units]' \
            + '\n'

    elif keywords[8] == 1:
        hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + date \
            + '\nOBJECT: ' + str(keywords[0]) \
            + '\nDATE-OBS: ' + str(keywords[1]) \
            + '\nRA: ' + str(keywords[2]) \
            + '\nDEC: ' + str(keywords[3]) \
            + '\nTELESCOPE: ' + str(keywords[4]) \
            + '\nINSTRUMENT: ' + str(keywords[5]) \
            + '\nEXPOSURE TIME: ' + str(keywords[6]) \
            + '\nID-SIM: ' + str(keywords[7]) \
            + '\nRV-FLAG: ' + str(keywords[8]) \
            + '\nRV: ' + str(keywords[9]) \
            + '\nRV-ERR: ' + str(keywords[10]) \
            + '\nRV-REF: ' + str(keywords[11]) \
            + '\nCOL1: WAVELENGH ARRAY [Angstrom]' \
            + '\nCOL2: NORMALISED FLUX ARRAY [arbitary units]' \
            + '\nCOL3: ERROR SPECTRUM [arbitary units]' \
            + '\n'

    return hascii


# Everything from here down is no longer used:


def get_header_data_aao(fib, path_final_vel, path_id_cat, comb_fits):
    """

    Return the necessary header data for AAO observations to be used in both
    ASCII and FITS files

    fib : int
        Integer denoting the AAO instrument fibre number

    path_final_vel : string
        Path to the appropriate "final_vel" file

    path_id_cat : string
        Path to the appropriate "master" file

    """

    # Load files needed for indexing    
    final_vel = ascii.read(path_final_vel)

    lmc_master = ascii.read(path_id_cat)

    vel_idx = np.where(fib == final_vel['FIBER'])[0][0]

    # comb_fits = fits.open(path_cfits)
    # Generate data to write
    date = get_date()
    sc = SkyCoord(final_vel['RA'][vel_idx], final_vel['DEC'][vel_idx],
                  unit=(u.hourangle, u.deg))
    obj_name_str = final_vel['NAME'][vel_idx].strip(' ')
    date_obs = comb_fits[2].header['ACTUTC']
    obj_ra = sc.ra.value
    obj_dec = sc.dec.value
    obj_tel = 'AAT'
    obj_ins = 'AAOmega'
    # obj_res = 'XXX'
    obj_rv = final_vel['VHEL'][vel_idx].strip(' ')
    texp = comb_fits[2].header['TWKODUR']
    print('{}.dat  {}   {}    {}    {}    {}    {}    {}'\
        .format(fib, obj_name_str, date_obs,
                obj_ra, obj_dec, obj_tel, obj_ins, texp))

    return date, obj_name_str, date_obs, obj_ra, obj_dec, obj_tel, obj_ins, obj_rv, texp


def get_date():
    """Return todays date in required ASTRO+ format""" 
    return str(datetime.datetime.today()).replace(" ", "T")[:-7]
