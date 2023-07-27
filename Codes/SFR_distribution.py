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

bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
bootesbestdf = pd.read_csv(bootesbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
bootesflux = pd.read_csv(bootesfluxfile, delimiter=',')

bootesdf = pd.merge(bootesdf, bootesbestdf, on=['name'], how='left')
bootesdf = pd.merge(bootesdf, bootesflux, on=['name'], how='left')

elaisfile = main + '/ELDF/File/en1_RGs_catalogue.csv'
elaisbest = main + '/ELDF/File/Crossmatch/full_en1_Best_crossmatch.csv'
elaisfluxfile = main + '/ELDF/File/en1_grg_flux.csv'

elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
elaisbestdf = pd.read_csv(elaisbest, delimiter=',', usecols=('name', 'AGN_final', 'Mass_cons','SFR_cons'))
elaisflux = pd.read_csv(elaisfluxfile, delimiter=',')

elaisdf = pd.merge(elaisdf, elaisbestdf, on=['name'], how='left')
elaisdf = pd.merge(elaisdf, elaisflux, on=['name'], how='left')

# df = df[df.lls < 0.7]

lhfile = main + '/LHLDF/File/lh_RGs_catalogue.csv'
lhbest = main + '/LHLDF/File/Crossmatch/full_lh_Best_crossmatch.csv'
lhfluxfile = main + '/LHLDF/File/lh_grg_flux.csv'

lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))
lhbestdf = pd.read_csv(lhbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
lhfluxdf = pd.read_csv(lhfluxfile, delimiter=',')

lhdf = pd.merge(lhdf, lhbestdf, on=['name'], how='left')
lhdf = pd.merge(lhdf, lhfluxdf, on=['name'], how='left')


df = pd.concat([bootesdf, elaisdf, lhdf])


df = df[~np.isnan(df.SFR_cons)]
df = df[df.SFR_cons > -10]

dfgrg = df[df.lls >= 0.7]
dfrg = df[df.lls < 0.7]

print(len(dfgrg.name), len(dfrg.name))
dfgrg = dfgrg[dfgrg['AGN_final']==1]
dfrg = dfrg[dfrg['AGN_final']==1]
print(len(dfgrg.name), len(dfrg.name))
exit()
print(ks_2samp(dfgrg.z, dfrg.z))

sns.histplot(dfrg.z, bins = 30, color = 'darkred', zorder=2, label = 'RGs')
sns.histplot(dfgrg.z, bins = 30, color = 'darkblue', zorder=3, label = 'GRGs', alpha = 0.5)
plt.show()
exit()

star_mass_grg = np.array(dfgrg.Mass_cons)
star_mass_grg = star_mass_grg[~np.isnan(star_mass_grg)]
star_mass_grg = star_mass_grg[star_mass_grg > 6]
sfr_grg = np.array(dfgrg.SFR_cons)
sfr_grg = sfr_grg[~np.isnan(sfr_grg)]
sfr_grg = sfr_grg[sfr_grg > -10]

star_mass_rg = np.array(dfrg.Mass_cons)
star_mass_rg = star_mass_rg[~np.isnan(star_mass_rg)]
star_mass_rg = star_mass_rg[star_mass_rg > 6]
sfr_rg = np.array(dfrg.SFR_cons)
sfr_rg = sfr_rg[~np.isnan(sfr_rg)]
sfr_rg = sfr_rg[sfr_rg > -10]


print(len(star_mass_rg), len(star_mass_grg))

# star_mass_grg = np.array(star_mass_grg)
# sfr_grg = np.array(sfr_grg)
#
# star_mass_rg = np.array(star_mass_rg)
# sfr_rg = np.array(sfr_rg)
# print(star_mass_grg)
binwidth = 0.35
xymax = max(np.max(np.abs(sfr_grg)), np.max(np.abs(sfr_rg)))
lim = (int(xymax/binwidth) + 1) * binwidth

print(len(sfr_rg))
print('sfr > 1:', (np.array(sfr_rg) > 1).sum())
# print(len(sfr_grg))

print(ks_2samp(sfr_rg, sfr_grg))
print(np.mean(sfr_rg), np.mean(sfr_grg))
print(np.std(sfr_rg), np.std(sfr_grg))

bins = np.arange(-lim, lim + binwidth, binwidth)
fig = plt.figure(figsize = (6, 6))
fig1, ax = plt.subplots()
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
# plt.plot(x, y, color = 'red', zorder = 4)
sns.histplot(sfr_rg, bins = bins, color = 'darkred', zorder=2, label = 'RGs')
sns.histplot(sfr_grg, bins = bins, color = 'darkblue', zorder=3, label = 'GRGs', alpha = 0.5)
plt.grid(color = 'white', zorder=0, linewidth = 1.5)
plt.xlabel(r'$Log(SFR[M_{\odot} yr^{-1}])$',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'N', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
ax.patch.set_facecolor('#ababab')
ax.patch.set_alpha(0.3)
plt.xlim(-3.0, 4.0)
plt.legend(fontsize=15)
plt.show()
