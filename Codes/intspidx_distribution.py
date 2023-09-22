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
from statistical_tests import StatisticalTests

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/File/lba-hba_int_spidx.txt'
rgfile = main + '/BLDF/File/bootes_RGs_catalogue.csv'


dfrg = pd.read_csv(rgfile, usecols=('name', 'lls'))
data = np.loadtxt(file, dtype={'names': ('col1', 'col2', 'col3'),
                                    'formats': ('U10', 'f8', 'f8')}, unpack=True)

df = pd.DataFrame({'name': data[0], 'spidx':data[1], 'spidxerr':data[2] })
df = pd.merge(df, dfrg[['name','lls']], on='name')



# number of RGs and GRGs
df_copy = df
df_copy = df_copy[df_copy.lls >= 0.7]
print("# of RGs:", len(df.name) - len(df_copy.name))
print("# of GRGs:", len(df_copy.name))

rgs_spidx = df.spidx[df.lls < 0.5]
grgs_spidx = df.spidx[df.lls >= 0.5]

test = StatisticalTests(rgs_spidx, grgs_spidx)

D,p_val = test.two_sample_ks_test()
print('Kolmogorv test:', D, p_val)
U,z,p_val = test.mann_whitney_u_test()
print('MW U test:', U,z, p_val)
exit()

# Create histograms using Seaborn
fig = plt.figure(figsize=(6, 6))  # Adjust the figure size as needed
fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

# Normalize histograms by dividing counts by the total number of elements

bins = np.histogram(np.concatenate((-rgs_spidx, -grgs_spidx)), bins=15)[1]
print(f'mean (rg):{np.mean(rgs_spidx)} (grg):{np.mean(grgs_spidx)}, median (rg):{np.median(rgs_spidx)} (grg):{np.median(grgs_spidx)}, 90perc. (rg):{np.percentile(rgs_spidx,5)}, (grg):{np.percentile(grgs_spidx,5)}')
sns.histplot(-grgs_spidx, bins=bins, kde=True, label=str(len(df_copy.name)) + ' GRGs', alpha=0.5, color='orange',stat='probability')
sns.histplot(-rgs_spidx, bins=bins, kde=True, label=str(len(df.name) - len(df_copy.name)) + ' RGs', alpha=0.5, color='blue', stat='probability')

# Add labels and title
plt.xlabel(r'$\alpha^{610}_{150}$', fontsize=15)
plt.ylabel('N', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend()
plt.show()


