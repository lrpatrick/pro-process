"""

PROMETEO functions to read input fits/ASCII files and create the standardised PROMETEO
FITS files


This script searches for keywords in the headers and if not avilable creates
these keywords. 


Author: LRP
Date: 13-02-2020
"""

import astropy.io.fits as fits
import astropy.io.ascii as ascii
import glob
import numpy as np
import os

import pfits

pid_start = int(input('[INFO] Please enter ID number for first spectrum as integer:\n')) or int(0)

# LMC
# tree = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/30_lmc_ccd2/'
# tree = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/2013_lmc_ccd2/'
# SMC
# tree = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/20_smc_ccd2/'
smc_dirs = ['20_smc_ccd2',  '21_smc_ccd2',  '28_smc_ccd2',  '29_smc_ccd2',  '30_smc_ccd2',  '6_smc_ccd2',  '7_smc_ccd2',  '8_smc_ccd2']

pid_fin = 0
for j, smc_dir in enumerate(smc_dirs):
    print('[INFO] Processing directory {}'.format(smc_dir))
    tree = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/' + smc_dir + '/'
    flist_aao = glob.glob(tree + '/[0-9][0-9][0-9].dat')

    # The three input files needed to create the headers:
    cfits = fits.open(tree + 'combined_frames.fits')

    # To be passed to pfits.write_fits() and pfits.write_ascii_header_aao()
    fv = tree + 'final_vel.cat'
    id_master = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/SMC'
    # id_master = '/media/lee/18FEB5E54AB6A1A6/data/PROMETEO/AAOmega/LMC'



    for i, faao in enumerate(flist_aao):
        print('[INFO] Processing file {}'.format(faao))
        pid_i = pid_start + i
        print('[INFO] PID {}'.format(pid_i))
        # Write ASCII
        faao_in = ascii.read(faao)
        fib = os.path.basename(faao)[:-4]
        try:
            aao_h1, obj_name_str = pfits.write_ascii_header_aao(fib, pid_i, final_vel=fv, id_cat=id_master, comb_fits=cfits)
            fname = obj_name_str + '_' + tree.split('/')[-2] + '_AAO_pro.dat'
            ascii_data = pfits.write_ascii(fname, faao_in['col1'], faao_in['col2'], head=aao_h1)
        except:
            print('[WARNING] ASCII file likely not in standard format')

        print('[INFO] Write output FITS file:')
        aao_outname_fits = obj_name_str + '_' + tree.split('/')[-2] + '_AAO_pro.fits'
        try:
            new_hdu = pfits.write_fits(cfits, ascii_data, pid_i, fib=fib, final_vel=fv, id_cat=id_master)
        except:
            print('[WARNING] FITS file likely not in standard format')


            # # Start with a fresh fits header based on the data
            # # Add important keywords:
            # new_hdu = pfits.write_fits(iacob_outname_fits, iacob1)
            # new_hdu = pfits.write_fits(iacob_outname_fits, iacob1, pid_i)


    pid_start = pid_i + 1
    print('[INFO] Final PID {}'.format(pid_i))