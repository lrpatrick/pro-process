"""

PROMETEO functions to read input data files and create the standardised
PROMETEO FITS and ASCII files


The project constsis of two main functions:

write_fits()
write_ascii()

These two functions are called from an input script and are used to write the
output files


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
    # 
    data_options = [x_option, y_option, yerr_option]
    data_keywords = [xname, x_ext, None, None, yname, y_ext]

    # Header:
    header_kws = ['OBJECT', 'DATE-OBS', 'RA', 'DEC',
                  'TELESCOP', 'INSTRUME', 'EXPTIME']
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
            i_match = [obj for obj in matches if 'I-' + kw == obj]
            if len(i_match) == 1:
                print('[INFO] IACOB keyword found')
                header_kws[i] = i_match[0]
                header_exts[i] = ext
            else:
                exact_match = [obj for obj in matches if kw == obj]
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

    fname = str(name) + '_' + str(date_obs) + '_' + str(tel) + '_' + str(ins)\
            + '_' + str(spec_start) + '_' + str(spec_end) + ext

    return fname


def write_fits_header(newfits, pid, keywords):
    """
    Create a fits header and instert it into the new fits fits
    
    pid should be an integer

    """
    # Generate date-time to put in header
    date = str(datetime.date.today())
    comment = 'ASTRO+ processed FITS file'
    # Put data into header
    newfits.header['DATE'] = date
    newfits.header['PROID'] = pid
    newfits.header['OBJECT'] = keywords[0]
    newfits.header['DATE-OBS'] = keywords[1]
    newfits.header['RA'] = keywords[2]
    newfits.header['DEC'] = keywords[3]
    newfits.header['TEL'] = keywords[4]
    newfits.header['INS'] = keywords[5]
    newfits.header['TEXP'] = keywords[6]
    newfits.header['COMMENT'] = comment

    return newfits


def write_fits(pid, data, keywords, oldfits=None):

    """
    Function to write the standard PROMETEO FITS files

    Arguments:

    pid : int
        ASTRO+ ID for star

    data : numpy.ndarray
        Array containing three columns:
            x : wavelength 
            y : flux 
            yerr : Uncertainity in flux

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
    
    oldfits : astropy.io.fits.hdu.hdulist.HDUList
        A fits file to be used to propagate header information

    Returns:
    prohdu : astropy.io.fits.hdu.hdulist.HDUList
        New FITS file
    """
    # Create clean FITS based on the data alone
    new_hdu = fits.PrimaryHDU(data)
    # Add header information to new fits
    new_hdu = write_fits_header(new_hdu, pid, keywords)
    # Get the filename to be written out
    fname = write_filename('.fits', keywords[0], keywords[1],
                                    keywords[4], keywords[5],
                                    int(round(data[:, 0][0])),
                                    int(round(data[:, 0][-1])))
    # Add additional information from FITS if appropriate
    if oldfits is not None:
        image_hdu = fits.ImageHDU(None, header=oldfits[0].header)
    else:
        # This is not a good solution!
        image_hdu = fits.ImageHDU(None)

    # Combine extensions into one file
    prohdu = fits.HDUList([new_hdu, image_hdu])

    # Write output FITS file
    prohdu.writeto('data/' + fname, output_verify='ignore')
    print('[INFO] FITS file written to {0}'.format(fname))
    return prohdu


def write_ascii(pid, data, keywords):
    """
    Write the standard PROMETEO ascii files
    data : numpy.ndarray
        A three column data array containing the following:
        x : numpy.ndarray
            Continually increasing wavelength axis for spectroscopic data
        y : numpy.ndarray
            Normalised spectroscopic data. Same size as x
        yerr : numpy.ndarray
            Uncertainity on y. Same size as x
    pid : int
        ASTRO+ ID for star

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



    """
    fname = write_filename('.dat', keywords[0], keywords[1],
                                    keywords[4], keywords[5],
                                    int(round(data[:, 0][0])),
                                    int(round(data[:, 0][-1])))
    head1 = write_ascii_header(pid, keywords)
    np.savetxt('data/' + fname, data, header=head1, fmt='%8.8f')
    print('[INFO] ASCII file written to {0}'.format(fname))
    return data


def write_ascii_header(pid, keywords):
    """
    Write ASCII standard header format

    """ 
    # Generate data to write
    date = str(datetime.date.today())
    pid = str(pid)

    hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + date \
        + '\nPROMETEO reference ID: ' + pid \
        + '\nOBJECT: ' + str(keywords[0]) \
        + '\nDATE-OBS: ' + str(keywords[1]) \
        + '\nRA: ' + str(keywords[2]) \
        + '\nDEC: ' + str(keywords[3]) \
        + '\nTELESCOPE: ' + str(keywords[4]) \
        + '\nINSTRUMENT: ' + str(keywords[5]) \
        + '\nEXPOSURE TIME: ' + str(keywords[6]) \
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
    date = str(datetime.date.today())
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


# def write_ascii_header_aao(fib, pid, final_vel, id_cat, comb_fits):
#     """
#     Write ASCII standard header format

#     """ 
#     # Generate header data
#     date, obj_name_str, date_obs, \
#           obj_ra, obj_dec, \
#           obj_tel, obj_ins, \
#           obj_rv, \
#           texp = get_header_data_aao(fib, final_vel, id_cat, comb_fits)

#     hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + str(date) \
#         + '\nPROMETEO reference ID: ' + str(pid) \
#         + '\nOBJECT: ' + str(obj_name_str) \
#         + '\nDATE-OBS: ' + str(date_obs) \
#         + '\nRA: ' + str(obj_ra) \
#         + '\nDEC: ' + str(obj_dec) \
#         + '\nTELESCOPE: ' + str(obj_tel) \
#         + '\nINSTRUMENT: ' + str(obj_ins) \
#         + '\nEXPOSURE TIME: ' + str(texp) \
#         + '\nHELIOCENTRIC RADIAL VELOCITY: ' + str(obj_rv) \
#         + '\nCOL1: WAVELENGH ARRAY [Angstrom]' \
#         + '\nCOL2: NORMALISED FLUX ARRAY [arbitary units]' \
#         + '\nCOL3: ERROR SPECTRUM [arbitary units]' \
#         + '\n'

#     return hascii, obj_name_str


# def write_fits_header_aao(newfits, pid, fib, final_vel, id_cat, comb_fits):
#     """
#     Create a fits header and instert it into the new fits fits
    
#     pid should be an integer

#     """
#     # Generate header data
#     date, obj_name_str, date_obs, \
#           obj_ra, obj_dec, \
#           obj_tel, obj_ins, \
#           obj_rv,\
#           texp = get_header_data_aao(fib, final_vel, id_cat, comb_fits)

#     # Put data into header
#     newfits.header['DATE'] = date
#     newfits.header['PROID'] = pid
#     newfits.header['OBJECT'] = obj_name_str
#     newfits.header['DATE-OBS'] = date_obs
#     # newfits.header['CRVAL1'] = crval1
#     # newfits.header['CDELT1'] = cdelt1

#     newfits.header['RA'] = obj_ra
#     newfits.header['DEC'] = obj_dec
#     newfits.header['TEL'] = obj_tel
#     newfits.header['INS'] = obj_ins
#     # newfits.header['RES'] = obj_res
#     newfits.header['VBAR'] = obj_rv
#     newfits.header['TEXP'] = texp

#     return newfits


# def get_header_data_harps(harps):
#     """

#     Return the necessary header data for AAO observations to be used in both
#     ASCII and FITS files

#     harps : astropy.fits.

#     """
#     date = str(datetime.date.today())
#     # pid = str()
#     obj_name_str =  harps[0].header['OBJECT']
#     date_obs = harps[0].header['DATE-OBS']
#     obj_ra = harps[0].header['RA']
#     obj_dec = harps[0].header['DEC']
#     obj_tel = harps[0].header['TELESCOP']
#     obj_ins = harps[0].header['INSTRUME']
#     texp = harps[0].header['EXPTIME']


#     return date, obj_name_str, date_obs, obj_ra, obj_dec, obj_tel, obj_ins, texp


# def write_fits_header_harps(oldfits, newfits, pid):
#     """
#     Create a fits header and instert it into the new fits fits
    
#     pid should be an integer

#     """
#     # Generate header data
#     date, obj_name_str, date_obs, \
#           obj_ra, obj_dec, \
#           obj_tel, obj_ins, \
#           texp = get_header_data_harps(oldfits)

#     # Put data into header
#     newfits.header['DATE'] = date
#     newfits.header['PROID'] = pid
#     newfits.header['OBJECT'] = obj_name_str
#     newfits.header['DATE-OBS'] = date_obs
#     newfits.header['RA'] = obj_ra
#     newfits.header['DEC'] = obj_dec
#     newfits.header['TEL'] = obj_tel
#     newfits.header['INS'] = obj_ins
#     newfits.header['TEXP'] = texp

#     return newfits


# def write_ascii_header_harps(fhdu, pid):
#     """
#     Write ASCII standard header format

#     """ 
#     # Generate header data
#     date, obj_name_str, date_obs, \
#           obj_ra, obj_dec, \
#           obj_tel, obj_ins, \
#           texp = get_header_data_harps(fhdu)


#     hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + str(date) \
#         + '\nPROMETEO reference ID: ' + str(pid) \
#         + '\nOBJECT: ' + str(obj_name_str) \
#         + '\nDATE-OBS: ' + str(date_obs) \
#         + '\nRA: ' + str(obj_ra) \
#         + '\nDEC: ' + str(obj_dec) \
#         + '\nTELESCOPE: ' + str(obj_tel) \
#         + '\nINSTRUMENT: ' + str(obj_ins) \
#         + '\nEXPOSURE TIME: ' + str(texp) \
#         + '\nCOL1: WAVELENGH ARRAY [Angstrom]' \
#         + '\nCOL2: NORMALISED FLUX ARRAY [arbitary units]' \
#         + '\nCOL3: ERROR SPECTRUM [arbitary units]' \
#         + '\n'

#     return hascii, obj_name_str


# def write_fits_harps(oldfits, data, pid):
#     """
#     Function to write the standard PROMETEO FITS files
    
#     This function current assumes that the data is the the zeroth extension

#     """
#     new_hdu = fits.PrimaryHDU(data)

#     new_hdu = write_fits_header_harps(oldfits, new_hdu, pid)

#     image_hdu = fits.ImageHDU(None, header=oldfits[0].header)

#     thdu = fits.HDUList([new_hdu, image_hdu])

#     fname = write_filename('.fits', new_hdu.header['OBJECT'],
#                                  new_hdu.header['DATE-OBS'],
#                                  new_hdu.header['TEL'],
#                                  new_hdu.header['INS'],
#                                  int(round(new_hdu.data[:, 0][0])),
#                                  int(round(new_hdu.data[:, 0][-1])))

#     thdu.writeto('data/' + fname, output_verify='ignore')
#     print('[INFO] FITS file written to {0}'.format(fname))
#     return thdu
