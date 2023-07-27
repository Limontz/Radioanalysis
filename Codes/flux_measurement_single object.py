from astropy.io import fits
from astropy.wcs import WCS
import astropy.utils as utils
import astropy
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5
import matplotlib.colors as colors
import pyregion
from function import SampleImage
from function import calculate_flux
from cutout2D import cutout
import astropy.units as u
import os, sys, argparse, logging
import pandas as pd

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
# fits_file = main + '/BLDF/File/bootes_radio_image.fits'
fits_file = main + '/BLDF/File/apt.fits'


name = 'J1424+3436'
#read the sample file
# sample_file = main + '/BLDF/File/bootes_GRGs_catalogue.csv'
# df = pd.read_csv(sample_file, delimiter=',')

# df.to_csv(main+'/ELDF/File/forvizier2.csv', columns = ['myra', 'mydec'], index=False)
# exit()



grgcutout, cdelt = cutout(fits_file, (float(216.1358))*u.deg, (float(34.6127))*u.deg, 45/60)
f1='hba.fits'
grgcutout.writeto(f1, overwrite=True)

regionfile = main + '/GMRT_proposal/DS9regions/' + name + '_rightbridge.reg'
# pyregion.open(regionfile)
flux, errors, upper = calculate_flux(f1, regionfile)

data = np.column_stack([ str(name), round(np.sum(flux) * 1e3,2) , round(np.sum(errors) *1e3,2)])
print(data)
# w_file = main + '/BLDF/File/bootes_grg_flux.txt'
# flux_file = open(w_file, "a+")
# np.savetxt(flux_file, data, fmt='%s')
# flux_file.close()
os.system('rm hba.fits')
