"""

PROMETEO functions to read input fits/ASCII files and create the standardised
PROMETEO FITS files


This script searches for keywords in the headers and if not avilable creates
these keywords. 

TODO:

Update filenames for ASCII files

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


def write_fits(oldfits, data, pid, fib=None, final_vel=None, id_cat=None):
    """
    Function to write the standard PROMETEO FITS files
    
    This function current assumes that the data is the the zeroth extension

    """
    new_hdu = fits.PrimaryHDU(data)
    if fib is None:
        new_hdu = write_fits_header(oldfits, new_hdu, pid)
    else:
        new_hdu = write_fits_header_aao(new_hdu, pid, fib, final_vel, id_cat, oldfits)

    fname = write_filename('.fits', new_hdu.header['OBJECT'],
                                 new_hdu.header['DATE-OBS'],
                                 new_hdu.header['TEL'],
                                 new_hdu.header['INS'],
                                 int(round(new_hdu.data[:, 0][0])),
                                 int(round(new_hdu.data[:, 0][-1])))

    image_hdu = fits.ImageHDU(None, header=oldfits[0].header)

    thdu = fits.HDUList([new_hdu, image_hdu])

    thdu.writeto('data/' + fname, output_verify='ignore')
    print('[INFO] FITS file written to {0}'.format(fname))
    return thdu


def write_fits_harps(oldfits, data, pid):
    """
    Function to write the standard PROMETEO FITS files
    
    This function current assumes that the data is the the zeroth extension

    """
    new_hdu = fits.PrimaryHDU(data)

    new_hdu = write_fits_header_harps(oldfits, new_hdu, pid)

    image_hdu = fits.ImageHDU(None, header=oldfits[0].header)

    thdu = fits.HDUList([new_hdu, image_hdu])

    fname = write_filename('.fits', new_hdu.header['OBJECT'],
                                 new_hdu.header['DATE-OBS'],
                                 new_hdu.header['TEL'],
                                 new_hdu.header['INS'],
                                 int(round(new_hdu.data[:, 0][0])),
                                 int(round(new_hdu.data[:, 0][-1])))

    thdu.writeto('data/' + fname, output_verify='ignore')
    print('[INFO] FITS file written to {0}'.format(fname))
    return thdu


def write_ascii(fname, x, y, yerr=0, head=''):
    """
    Write the standard PROMETEO ascii files
    
    fname : string
        Name for the ASCII file to be written. Should end in '.dat'

    x : numpy.ndarray
        Continually increasing wavelength axis for spectroscopic data

    y : numpy.ndarray
        Normalised spectroscopic data. Same size as x

    head : string
        header information to be written to file in standard format 

    """
    # print(yerr.shape[0] - y.shape[0])
    # if yerr.shape[0] != y.shape[0]:
    yerr = np.zeros_like(x)
        # print('Yerr All zeros')
    out_data = np.column_stack((x, y, yerr))
    np.savetxt('data/' + fname, out_data, header=head, fmt='%8.8f')
    print('[INFO] ASCII file written to {0}'.format(fname))
    return out_data



def write_fits_header(oldfits, newfits, pid):
    """
    Create a fits header and instert it into the new fits fits
    
    pid should be an integer

    """
    # Generate header data
    date = str(datetime.date.today())
    # pid = str()
    obj_name_str =  oldfits[0].header['OBJECT']
    date_obs = oldfits[0].header['DATE-OBS']
    obj_ra = oldfits[0].header['I-RA']
    obj_dec = oldfits[0].header['I-DEC']
    crval1 = oldfits[0].header['CRVAL1']
    cdelt1 = oldfits[0].header['CDELT1']
    obj_tel = oldfits[0].header['TELESCOP']
    obj_ins = oldfits[0].header['INSTRUME']
    # obj_res = oldfits[0].header['I-RESOL']
    obj_rv = oldfits[0].header['I-VBAR']
    texp = oldfits[0].header['I-TEXP']
    # Put data into header
    newfits.header['DATE'] = date
    newfits.header['PROID'] = pid
    newfits.header['OBJECT'] = obj_name_str
    newfits.header['DATE-OBS'] = date_obs
    # newfits.header['CRVAL1'] = crval1
    # newfits.header['CDELT1'] = cdelt1

    newfits.header['RA'] = obj_ra
    newfits.header['DEC'] = obj_dec
    newfits.header['TEL'] = obj_tel
    newfits.header['INS'] = obj_ins
    # newfits.header['RES'] = obj_res
    newfits.header['VBAR'] = obj_rv
    newfits.header['TEXP'] = texp

    return newfits


def write_ascii_header(fhdu, pid):
    """
    Write ASCII standard header format

    """ 
    # Generate data to write
    date = str(datetime.date.today())
    pid = str(pid)
    obj_name_str =  fhdu[0].header['OBJECT']
    date_obs = fhdu[0].header['DATE-OBS']
    obj_ra_str = str(fhdu[0].header['RA'])
    obj_dec_str = str(fhdu[0].header['DEC'])
    obj_tel = str(fhdu[0].header['TELESCOP'])
    obj_ins = str(fhdu[0].header['INSTRUME'])
    # obj_res = str(fhdu[0].header['I-RESOL'])
    obj_rv = str(fhdu[0].header['I-VBAR'])
    texp = str(fhdu[0].header['I-TEXP'])

    hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + date \
        + '\nPROMETEO reference ID: ' + pid \
        + '\nOBJECT: ' + obj_name_str \
        + '\nDATE-OBS: ' + date_obs \
        + '\nRA: ' + obj_ra_str \
        + '\nDEC: ' + obj_dec_str \
        + '\nTELESCOPE: ' + obj_tel \
        + '\nINSTRUMENT: ' + obj_ins \
        + '\nEXPOSURE TIME: ' + str(texp) \
        + '\nHELIOCENTRIC RADIAL VELOCITY: ' + obj_rv \
        + '\nCOL1: WAVELENGH ARRAY [Angstrom]' \
        + '\nCOL2: NORMALISED FLUX ARRAY [arbitary units]' \
        + '\nCOL3: ERROR SPECTRUM [arbitary units]' \
        + '\n'
        # + '\nSPECTRAL RESOLVING POWER: ' + obj_res \

    return hascii


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

    return date, obj_name_str, date_obs, obj_ra, obj_dec, obj_tel, obj_ins, obj_rv, texp


def write_ascii_header_aao(fib, pid, final_vel, id_cat, comb_fits):
    """
    Write ASCII standard header format

    """ 
    # Generate header data
    date, obj_name_str, date_obs, \
          obj_ra, obj_dec, \
          obj_tel, obj_ins, \
          obj_rv, \
          texp = get_header_data_aao(fib, final_vel, id_cat, comb_fits)

    hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + str(date) \
        + '\nPROMETEO reference ID: ' + str(pid) \
        + '\nOBJECT: ' + str(obj_name_str) \
        + '\nDATE-OBS: ' + str(date_obs) \
        + '\nRA: ' + str(obj_ra) \
        + '\nDEC: ' + str(obj_dec) \
        + '\nTELESCOPE: ' + str(obj_tel) \
        + '\nINSTRUMENT: ' + str(obj_ins) \
        + '\nEXPOSURE TIME: ' + str(texp) \
        + '\nHELIOCENTRIC RADIAL VELOCITY: ' + str(obj_rv) \
        + '\nCOL1: WAVELENGH ARRAY [Angstrom]' \
        + '\nCOL2: NORMALISED FLUX ARRAY [arbitary units]' \
        + '\nCOL3: ERROR SPECTRUM [arbitary units]' \
        + '\n'

    return hascii, obj_name_str


def write_fits_header_aao(newfits, pid, fib, final_vel, id_cat, comb_fits):
    """
    Create a fits header and instert it into the new fits fits
    
    pid should be an integer

    """
    # Generate header data
    date, obj_name_str, date_obs, \
          obj_ra, obj_dec, \
          obj_tel, obj_ins, \
          obj_rv,\
          texp = get_header_data_aao(fib, final_vel, id_cat, comb_fits)

    # Put data into header
    newfits.header['DATE'] = date
    newfits.header['PROID'] = pid
    newfits.header['OBJECT'] = obj_name_str
    newfits.header['DATE-OBS'] = date_obs
    # newfits.header['CRVAL1'] = crval1
    # newfits.header['CDELT1'] = cdelt1

    newfits.header['RA'] = obj_ra
    newfits.header['DEC'] = obj_dec
    newfits.header['TEL'] = obj_tel
    newfits.header['INS'] = obj_ins
    # newfits.header['RES'] = obj_res
    newfits.header['VBAR'] = obj_rv
    newfits.header['TEXP'] = texp

    return newfits


def get_header_data_harps(harps):
    """

    Return the necessary header data for AAO observations to be used in both
    ASCII and FITS files

    harps : astropy.fits.

    """
    date = str(datetime.date.today())
    # pid = str()
    obj_name_str =  harps[0].header['OBJECT']
    date_obs = harps[0].header['DATE-OBS']
    obj_ra = harps[0].header['RA']
    obj_dec = harps[0].header['DEC']
    obj_tel = harps[0].header['TELESCOP']
    obj_ins = harps[0].header['INSTRUME']
    texp = harps[0].header['EXPTIME']


    return date, obj_name_str, date_obs, obj_ra, obj_dec, obj_tel, obj_ins, texp


def write_fits_header_harps(oldfits, newfits, pid):
    """
    Create a fits header and instert it into the new fits fits
    
    pid should be an integer

    """
    # Generate header data
    date, obj_name_str, date_obs, \
          obj_ra, obj_dec, \
          obj_tel, obj_ins, \
          texp = get_header_data_harps(oldfits)

    # Put data into header
    newfits.header['DATE'] = date
    newfits.header['PROID'] = pid
    newfits.header['OBJECT'] = obj_name_str
    newfits.header['DATE-OBS'] = date_obs
    newfits.header['RA'] = obj_ra
    newfits.header['DEC'] = obj_dec
    newfits.header['TEL'] = obj_tel
    newfits.header['INS'] = obj_ins
    newfits.header['TEXP'] = texp

    return newfits


def write_ascii_header_harps(fhdu, pid):
    """
    Write ASCII standard header format

    """ 
    # Generate header data
    date, obj_name_str, date_obs, \
          obj_ra, obj_dec, \
          obj_tel, obj_ins, \
          texp = get_header_data_harps(fhdu)


    hascii = 'PROMETEO processed ASCII file\nProcessing date: ' + str(date) \
        + '\nPROMETEO reference ID: ' + str(pid) \
        + '\nOBJECT: ' + str(obj_name_str) \
        + '\nDATE-OBS: ' + str(date_obs) \
        + '\nRA: ' + str(obj_ra) \
        + '\nDEC: ' + str(obj_dec) \
        + '\nTELESCOPE: ' + str(obj_tel) \
        + '\nINSTRUMENT: ' + str(obj_ins) \
        + '\nEXPOSURE TIME: ' + str(texp) \
        + '\nCOL1: WAVELENGH ARRAY [Angstrom]' \
        + '\nCOL2: NORMALISED FLUX ARRAY [arbitary units]' \
        + '\nCOL3: ERROR SPECTRUM [arbitary units]' \
        + '\n'

    return hascii, obj_name_str
