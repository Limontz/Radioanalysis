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
from scipy.stats import ks_2samp
#from astropy.cosmology import Planck13 as cosmo
from astropy import wcs
from astropy.io import fits
from astropy import constants as c
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
cosmo =  FlatLambdaCDM(H0=70, Om0=0.3)
import scipy.stats as sts

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/File/gmrt-hba_int_spidx.txt'
rgfile = main + '/ELDF/File/en1_RGs_catalogue.csv'


dfrg = pd.read_csv(rgfile, usecols=('name', 'lls', 'z'))
data = np.loadtxt(file, dtype={'names': ('col1', 'col2', 'col3'),
                                    'formats': ('U10', 'f8', 'f8')}, unpack=True)

df = pd.DataFrame({'name': data[0], 'spidx':data[1], 'spidxerr':data[2] })
df = pd.merge(df, dfrg[['name','lls', 'z']], on='name')


# number of RGs and GRGs
df_copy = df
df_copy = df_copy[df_copy.lls >= 0.7]
print("# of RGs:", len(df.name) - len(df_copy.name))
print("# of GRGs:", len(df_copy.name))

print(sts.pearsonr(df.lls,-df.spidx))

# Create histograms using Seaborn
fig = plt.figure(figsize=(7,6))  # Adjust the figure size as needed
fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
plt.errorbar(df.lls, -df.spidx, yerr=df.spidxerr, fmt='', capsize=1, elinewidth=0.2, color='gray', linestyle="")
plt.scatter(df.lls, -df.spidx, marker = ".", s= 100, c=df.z, cmap='rainbow', vmax=1.5)
#colorbar
colorbar = plt.colorbar()
colorbar.set_label('z', fontsize=15)
colorbar.ax.tick_params(labelsize=15)

xticks = [0.2, 0.6, 1.0, 2.0]

# Add labels and title
plt.ylabel(r'$\alpha^{610}_{150}$', fontsize=15)
plt.xlabel('LLS', fontsize=15)
plt.xscale('log')
plt.xticks(xticks, xticks, fontsize=15)
plt.yticks(fontsize=15)

plt.show()





