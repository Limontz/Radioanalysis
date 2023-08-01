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
from function import flatten

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/LHLDF/File/apt_LH_rgs.csv'
grgfile = main + '/LHLDF/File/LH_RGs_catalogue.csv'
hbafile = main + '/LHLDF/File/lockman_radio_image.fits'
aptfile = main + '/LHLDF/File/LH_LofarApertif.fits'
# lbafile = main + '/BLDF/File/Bootes_LoLLS.fits'
# gmrtfile = main + '/ELDF/File/GMRT_EN1_dr1.fits'

freq_name = 'apt'

dfgrg = pd.read_csv(grgfile, delimiter=",", usecols=('name', 'ra', 'dec', 'las', 'z', 'lls'))
df = pd.read_csv(file)

df = pd.merge(df['name'], dfgrg[['name','ra', 'dec', 'las']], on='name')

for i in range(len(df.name)):
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 4D
    # HBA MAPS
    hdu4D = fits.open(hbafile)



    header, data = flatten(hbafile)

    new_fh = fits.PrimaryHDU(data=data, header=header)
    path=main + '/LHLDF/File/spidxmap/hba/'+str(df.name[i])+'_hba.fits'
    new_fitsfile = str(df.name[i])+'_hba.fits'
    new_fh.writeto(new_fitsfile, overwrite=True)
    file = cutout(new_fitsfile, df.ra[i]*u.deg, df.dec[i]*u.deg, 6*df.las[i]/60, path)
    os.remove(new_fitsfile)
    # print('----------------------------------------------------------------------------------------------')



    # new_header = mywcs.to_header()
    header, data = flatten(aptfile)
    # header['BMAJ'] = 6/3600
    # header['BMIN'] = 5/3600
    # header['BPA'] = 45
    new_fh = fits.PrimaryHDU(data=data, header=header)
    path= main + '/LHLDF/File/spidxmap/apt/' +str(df.name[i])+'_'+freq_name+'.fits'
    new_fitsfile = str(df.name[i])+freq_name
    new_fh.writeto(new_fitsfile, overwrite=True)
    
    file = cutout(new_fitsfile, df.ra[i]*u.deg, df.dec[i]*u.deg, 6*df.las[i]/60, path)
    os.remove(new_fitsfile)
    # print('----------------------------------------------------------------------------------------------')


    # header, data = flatten(aptfile)

    # new_fh = fits.PrimaryHDU(data=data, header=header)
    # new_fitsfile='2D-'+ 'apt.fits'
    # new_fh.writeto(new_fitsfile, overwrite=True)

    # file = cutout(new_fitsfile, dfgrg.ra[i]*u.deg, dfgrg.dec[i]*u.deg, 6*dfgrg.las[i]/60, dfgrg.name[i], 'BLDF', 'apt')
    # os.remove(new_fitsfile)