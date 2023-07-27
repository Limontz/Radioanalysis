import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

import numpy as  np
import pandas as pd
import matplotlib.font_manager

from uncertainties import unumpy
from scipy import optimize

import aplpy
import pyregion
from astropy.coordinates import SkyCoord
from regions import PixCoord
from regions import read_ds9
#from astropy.cosmology import Planck13 as cosmo
from astropy import wcs
from astropy.io import fits
from astropy import constants as c
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cutout2D_creation import cutout

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
grgfile = main + '/BLDF/File/bootes_GRGs_catalogue.csv'
hbafile = main + '/BLDF/File/bootes_radio_image.fits'
aptfile = main + '/BLDF/File/apt.fits'
lbafile = main + '/BLDF/File/Bootes_LoLLS.fits'

dfgrg = pd.read_csv(grgfile, delimiter=",", usecols=('name', 'ra', 'dec', 'las', 'z', 'lls'))
i = 14
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 4D
# HBA MAPS
hdu4D = fits.open(hbafile)
# hdu4D.info()

from function import flatten
header, data = flatten(hbafile)
# data = hdu4D[0].data.squeeze() # drops the size-1 axes
# data = data * 1.e3
# header = hdu4D[0].header
# # print(header)
# mywcs = wcs.WCS(header).celestial

# new_header = mywcs.to_header()

new_fh = fits.PrimaryHDU(data=data, header=header)
new_fitsfile='2D-'+ 'hba.fits'
new_fh.writeto(new_fitsfile, overwrite=True)

file = cutout(new_fitsfile, dfgrg.ra[i]*u.deg, dfgrg.dec[i]*u.deg, 6*dfgrg.las[i]/60, dfgrg.name[i], 'BLDF', 'hba')
os.remove(new_fitsfile)
print('----------------------------------------------------------------------------------------------')

# hdu4D = fits.open(lbafile)
# hdu4D.info()


# data = hdu4D[0].data.squeeze() # drops the size-1 axes
# data = data * 1.e3
# header = hdu4D[0].header
# #print(header)
# mywcs = wcs.WCS(header).celestial

# new_header = mywcs.to_header()
header, data = flatten(lbafile)

new_fh = fits.PrimaryHDU(data=data, header=header)
new_fitsfile='2D-'+ 'lba.fits'
new_fh.writeto(new_fitsfile, overwrite=True)

file = cutout(new_fitsfile, dfgrg.ra[i]*u.deg, dfgrg.dec[i]*u.deg, 6*dfgrg.las[i]/60, dfgrg.name[i], 'BLDF', 'lba')
os.remove(new_fitsfile)
print('----------------------------------------------------------------------------------------------')

# hdu4D = fits.open(aptfile)
# hdu4D.info()


# data = hdu4D[0].data.squeeze() # drops the size-1 axes
# data = data * 1.e3
# header = hdu4D[0].header
# #print(header)
# mywcs = wcs.WCS(header).celestial

# new_header = mywcs.to_header()
header, data = flatten(aptfile)

new_fh = fits.PrimaryHDU(data=data, header=header)
new_fitsfile='2D-'+ 'apt.fits'
new_fh.writeto(new_fitsfile, overwrite=True)

file = cutout(new_fitsfile, dfgrg.ra[i]*u.deg, dfgrg.dec[i]*u.deg, 6*dfgrg.las[i]/60, dfgrg.name[i], 'BLDF', 'apt')
os.remove(new_fitsfile)