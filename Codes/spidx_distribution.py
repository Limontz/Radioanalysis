import numpy as  np
import pandas as pd
import matplotlib.font_manager
import matplotlib.pyplot as plt
from uncertainties import unumpy
from scipy import optimize
from astropy.io import fits
import seaborn as sns
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
cosmo =  FlatLambdaCDM(H0=70, Om0=0.3)

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/ELDF/File/GMRT_en1_RGs.csv'
rgfile = main + '/ELDF/File/en1_RGs_catalogue.csv'

df = pd.read_csv(file)
dfrg = pd.read_csv(rgfile)

df = pd.merge(df['name'], dfrg[['name','lls']], on='name')
df_copy = df

# number of RGs and GRGs
df_copy = df_copy[df_copy.lls >= 0.7]
print("# of RGs:", len(df.name) - len(df_copy.name))
print("# of GRGs:", len(df_copy.name))


rgs_spidx = np.array([])
grgs_spidx = np.array([])
for i in range(len(df.name)):

    spidx_file = main + '/ELDF/File/spidxmap/output/gmrt-hba/' + str(df.name[i]) + '_spidx.fits'
    spidx_err_file = main + '/ELDF/File/spidxmap/output/gmrt-hba/' + str(df.name[i]) + '_spidx-err.fits'

    # Open the FITS file
    with fits.open(spidx_file) as hdul:
         # Assuming the data is in the first HDU (header data unit)
         data = hdul[0].data
    with fits.open(spidx_err_file) as hdul:
         # Assuming the data is in the first HDU (header data unit)
         data_err = hdul[0].data

    # Remove NaN values and create a new array without them
    data_cleaned = data[~np.isnan(data)]
    data_err_cleaned = data_err[~np.isnan(data)]

    data = data_cleaned[data_err_cleaned <= 0.15] #0.15 for GMRT data
    data_err = data_err_cleaned[data_err_cleaned <= 0.15]


    if (df.lls[i] < 0.7): 
        rgs_spidx = np.concatenate((rgs_spidx, data))
    else:
        grgs_spidx = np.concatenate((grgs_spidx, data))


# Create histograms using Seaborn
fig = plt.figure(figsize=(6, 6))  # Adjust the figure size as needed
fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

# Normalize histograms by dividing counts by the total number of elements

bins = np.histogram(np.concatenate((rgs_spidx, grgs_spidx)), bins=20)[1]

  
sns.histplot(-grgs_spidx, bins=bins, kde=True, label='GRGs', alpha=0.5, color='orange', stat="probability")
sns.histplot(-rgs_spidx, bins=bins, kde=True, label='RGs', alpha=0.5, color='blue', stat="probability")

# Add labels and title
plt.xlabel(r'$\alpha$', fontsize=15)
plt.ylabel('Normalised frequency', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend()
plt.show()
