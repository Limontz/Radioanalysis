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
#from astropy.cosmology import Planck13 as cosmo
from astropy.cosmology import LambdaCDM
from astropy import units as u
import seaborn as sns
import matplotlib.image as mpimg
from PIL import Image
import pandas as pd
from scipy.stats import ks_2samp

def get_data(fname,colname):
    data=fits.open(fname)
    data=data[1].data
    return data[colname]


main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
bootesfile = main + '/BLDF/File/bootes_RGs_catalogue.csv'
bootesbest = main + '/BLDF/File/Crossmatch/full_Bootes_Best_crossmatch.csv'
bootesfluxfile = main + '/BLDF/File/bootes_GRG_flux.csv'

<<<<<<< HEAD
bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'ztype', 'lls'))
=======
bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
>>>>>>> a06345a (Initial commit)
bootesbestdf = pd.read_csv(bootesbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
bootesflux = pd.read_csv(bootesfluxfile, delimiter=',')

bootesdf = pd.merge(bootesdf, bootesbestdf, on=['name'], how='left')
bootesdf = pd.merge(bootesdf, bootesflux, on=['name'], how='left')

elaisfile = main + '/ELDF/File/en1_RGs_catalogue.csv'
elaisbest = main + '/ELDF/File/Crossmatch/full_en1_Best_crossmatch.csv'
elaisfluxfile = main + '/ELDF/File/en1_grg_flux.csv'

<<<<<<< HEAD
elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'ztype', 'lls'))
=======
elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
>>>>>>> a06345a (Initial commit)
elaisbestdf = pd.read_csv(elaisbest, delimiter=',', usecols=('name', 'AGN_final', 'Mass_cons','SFR_cons'))
elaisflux = pd.read_csv(elaisfluxfile, delimiter=',')

elaisdf = pd.merge(elaisdf, elaisbestdf, on=['name'], how='left')
elaisdf = pd.merge(elaisdf, elaisflux, on=['name'], how='left')

# df = df[df.lls < 0.7]

lhfile = main + '/LHLDF/File/lh_RGs_catalogue.csv'
lhbest = main + '/LHLDF/File/Crossmatch/full_lh_Best_crossmatch.csv'
lhfluxfile = main + '/LHLDF/File/lh_grg_flux.csv'

<<<<<<< HEAD
lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'ztype', 'lls'))
=======
lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))
>>>>>>> a06345a (Initial commit)
lhbestdf = pd.read_csv(lhbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
lhfluxdf = pd.read_csv(lhfluxfile, delimiter=',')

lhdf = pd.merge(lhdf, lhbestdf, on=['name'], how='left')
lhdf = pd.merge(lhdf, lhfluxdf, on=['name'], how='left')


df = pd.concat([bootesdf, elaisdf, lhdf])

<<<<<<< HEAD
df = df[df.ztype=='s']

=======
>>>>>>> a06345a (Initial commit)
dfgrg = df[df.lls >= 0.7]
dfrg = df[df.lls < 0.7]


star_mass_grg = np.array(dfgrg.Mass_cons)
star_mass_grg = star_mass_grg[~np.isnan(star_mass_grg)]
sfr_grg = np.array(dfgrg.SFR_cons)
sfr_grg = sfr_grg[~np.isnan(sfr_grg)]

star_mass_rg = np.array(dfrg.Mass_cons)
star_mass_rg = star_mass_rg[~np.isnan(star_mass_rg)]
sfr_rg = np.array(dfrg.SFR_cons)
sfr_rg = sfr_rg[~np.isnan(sfr_rg)]


print(len(star_mass_rg), len(star_mass_grg))

# star_mass_grg = np.array(star_mass_grg)
# sfr_grg = np.array(sfr_grg)
#
# star_mass_rg = np.array(star_mass_rg)
# sfr_rg = np.array(sfr_rg)
# print(star_mass_grg)
binwidth = 0.20
xymax = max(np.max(np.abs(star_mass_grg)), np.max(np.abs(star_mass_rg)))
lim = (int(xymax/binwidth) + 1) * binwidth

print(len(sfr_grg))
print('sfr > 1:', (np.array(sfr_grg) > 1).sum())
print(len(sfr_grg))

print(ks_2samp(star_mass_rg, star_mass_grg))

bins = np.arange(-lim, lim + binwidth, binwidth)
fig = plt.figure(figsize = (6, 6))
fig1, ax = plt.subplots()
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
# plt.plot(x, y, color = 'red', zorder = 4)
sns.histplot(star_mass_rg, bins = bins,edgecolor='black', linewidth=1, color = 'darkred', zorder=2, label = 'RGs')
sns.histplot(star_mass_grg, bins = bins, edgecolor='black', linewidth=1,  color = 'darkblue', zorder=3, label = 'GRGs', alpha = 0.5)
plt.grid(color = 'white', zorder=0, linewidth = 1.5)
plt.xlabel(r'$Log(M_*/M_{\odot})$',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'N', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
ax.patch.set_facecolor('#ababab')
ax.patch.set_alpha(0.3)
plt.xlim(9.0, 13.0)
plt.legend(fontsize=15)
plt.show()
