import numpy as np
import astropy
import pandas as pd
import datashader as ds
import colorcet as cc
#from astropy.cosmology import Planck13 as cosmo
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
import matplotlib.pyplot as plt
from lib_linearfit import linear_fit, linear_fit_bootstrap
from astropy.coordinates import SkyCoord
from astropy import units as u
import seaborn as sns
from matplotlib.ticker import ScalarFormatter, NullFormatter
from mpl_toolkits.mplot3d import Axes3D
import time
from scipy.stats import ks_2samp

Mpc_to_m = 3.08e22
Jy_to_W = 1.e-26
alpha = 0.7

def search_around(coord1, coord2, z1, z2, seplimit):
    c1 = SkyCoord(ra=coord1[0]*u.deg, dec=coord1[1]*u.deg, distance=cosmo.comoving_distance(z1), frame='icrs')
    # c2 = SkyCoord(ra=datadesi['ra'].head(len)*u.deg, dec=datadesi['dec'].head(len)*u.deg, distance=cosmo.comoving_distance(datadesi['z_phot_mean'].head(len)), frame='icrs')
    c2 = SkyCoord(ra=coord2[0]*u.deg, dec=coord2[1]*u.deg, distance=cosmo.comoving_distance(z2), frame='icrs')
    idx2, idx1, sep2d, sep3d = SkyCoord.search_around_3d(c1, c2, seplimit*u.Mpc)
    return(idx2, idx1, sep2d, sep3d)

#******************* file reading and data management ****************************************

# Bootes field
main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
file = main + 'BLDF/File/bootes_GRGs_catalogue.csv'
desi_file = main + 'BLDF/File/DESI_Bootes-photo-z.csv'
bootesbest = main + '/BLDF/File/Crossmatch/full_Bootes_Best_crossmatch.csv'
bootesfluxfile = main + '/BLDF/File/bootes_GRG_flux.csv'

bootesdf = pd.read_csv(file, delimiter=None, usecols=('name', 'ra', 'dec', 'z', 'zstd', 'lls'))
bootesdesi = pd.read_csv(desi_file, delimiter=None, usecols=('ra', 'dec', 'dered_mag_r', 'z_phot_mean', 'z_spec', 'z_phot_std'))
bootesbestdf = pd.read_csv(bootesbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
bootesflux = pd.read_csv(bootesfluxfile, delimiter=',')

bootesdf = pd.merge(bootesdf, bootesflux, on=['name'], how='left')
bootesdf['power'] = 4. * np.pi * (1 + bootesdf.z)**(alpha - 1) * (bootesdf.flux/1e3) *  ((cosmo.luminosity_distance(bootesdf.z) * Mpc_to_m / u.Mpc)**2) * Jy_to_W   #units W


bootesdesi['z_phot_mean'] = np.where(bootesdesi['z_spec'] > 0, bootesdesi['z_spec'], bootesdesi['z_phot_mean'])
bootesdesi = bootesdesi[bootesdesi['z_phot_mean']>0]


# dfgrg_copy = dfgrg
# datadesi_copy = datadesi

# Elais field
elaisfile = main + '/ELDF/File/en1_GRGs_catalogue.csv'
desifile = main + '/ELDF/FILE/Elais-DESI-photo-z.csv'
elaisbest = main + '/ELDF/File/Crossmatch/full_en1_Best_crossmatch.csv'
elaisfluxfile = main + '/ELDF/File/en1_grg_flux.csv'

elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'ra', 'dec', 'z', 'zstd', 'lls'))
elaisdesi = pd.read_csv(desifile, delimiter=',', usecols=('ra', 'dec', 'dered_mag_r', 'z_phot_mean', 'z_spec', 'z_phot_std'))
elaisbestdf = pd.read_csv(elaisbest, delimiter=',', usecols=('name', 'AGN_final', 'Mass_cons','SFR_cons'))
elaisflux = pd.read_csv(elaisfluxfile, delimiter=',')

elaisdf = pd.merge(elaisdf, elaisflux, on=['name'], how='left')
elaisdf['power'] = 4. * np.pi * (1 + elaisdf.z)**(alpha - 1) * (elaisdf.flux/1e3) *  ((cosmo.luminosity_distance(elaisdf.z) * Mpc_to_m / u.Mpc)**2) * Jy_to_W   #units W

elaisdesi['z_phot_mean'] = np.where(elaisdesi['z_spec'] > 0, elaisdesi['z_spec'], elaisdesi['z_phot_mean'])
elaisdesi = elaisdesi[elaisdesi['z_phot_mean']>0]

# Lockman field 

lhfile = main + '/LHLDF/File/lh_GRGs_catalogue.csv'
desifile = main + '/LHLDF/FILE/DESI-LH_photo-z.csv'
lhbest = main + '/LHLDF/File/Crossmatch/full_lh_Best_crossmatch.csv'
lhfluxfile = main + '/LHLDF/File/lh_grg_flux.csv'

lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'ra', 'dec', 'z', 'zstd', 'lls'))
lhdesi = pd.read_csv(desifile, delimiter=',', usecols=('ra', 'dec', 'dered_mag_r', 'z_phot_mean', 'z_spec', 'z_phot_std'))
lhbestdf = pd.read_csv(lhbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
lhfluxdf = pd.read_csv(lhfluxfile, delimiter=',')

lhdf = pd.merge(lhdf, lhfluxdf, on=['name'], how='left')
lhdf['power'] = 4. * np.pi * (1 + lhdf.z)**(alpha - 1) * (lhdf.flux/1e3) *  ((cosmo.luminosity_distance(lhdf.z) * Mpc_to_m / u.Mpc)**2) * Jy_to_W   #units W


lhdesi['z_phot_mean'] = np.where(lhdesi['z_spec'] > 0, lhdesi['z_spec'], lhdesi['z_phot_mean'])
lhdesi = lhdesi[lhdesi['z_phot_mean']>0]


depthg = 24.0
depthr = 23.4
depthz = 22.5

zmax = 0.7 #np.max(dfgrg['z'])  #
elaisdf = elaisdf[elaisdf['z'] <= zmax]
bootesdf = bootesdf[bootesdf['z'] <= zmax]
lhdf = lhdf[lhdf['z'] <= zmax]


depthg = 24.0
depthr = 23.4
depthz = 22.5

mdesi_bootes = depthr - 5 * np.log10(cosmo.luminosity_distance(zmax).to(u.pc)/u.pc) + 5 * np.log10(cosmo.luminosity_distance(bootesdesi['z_phot_mean']).to(u.pc)/u.pc) #+Kcorrection?
mdesi_elais = depthr - 5 * np.log10(cosmo.luminosity_distance(zmax).to(u.pc)/u.pc) + 5 * np.log10(cosmo.luminosity_distance(elaisdesi['z_phot_mean']).to(u.pc)/u.pc) #+Kcorrection?
mdesi_lockman = depthr - 5 * np.log10(cosmo.luminosity_distance(zmax).to(u.pc)/u.pc) + 5 * np.log10(cosmo.luminosity_distance(lhdesi['z_phot_mean']).to(u.pc)/u.pc) #+Kcorrection?

elaisdesi = elaisdesi[elaisdesi['dered_mag_r'] <= mdesi_elais]
lhdesi = lhdesi[lhdesi['dered_mag_r'] <= mdesi_lockman]
bootesdesi = bootesdesi[bootesdesi['dered_mag_r'] <= mdesi_bootes]

#*********************************************************************************************************************************************************
#******************* number of galaxies around the grgs **************************************************************************************************
#*********************************************************************************************************************************************************

maxdistance = 10 #Mpc unit

# Bootes field
idxdesi, idxgrg, sep2d, sep3d = search_around((bootesdf['ra'], bootesdf['dec']), (bootesdesi['ra'], bootesdesi['dec']), bootesdf['z'], bootesdesi['z_phot_mean'], maxdistance)

bootes_dftotal = pd.DataFrame({'name':np.array(bootesdf['name'])[idxgrg], 'ra':np.array(bootesdf['ra'])[idxgrg], 'dec':np.array(bootesdf['dec'])[idxgrg], 'z':np.array(bootesdf['z'])[idxgrg], 'radesi':bootesdesi.iloc[idxdesi, 0],\
                        'decdesi':bootesdesi.iloc[idxdesi, 1], 'sep3d': sep3d})

table = bootes_dftotal['name'].value_counts()
bootes_df = pd.DataFrame({'name':table.index, 'ngalaxies':table.values})
bootes_df = pd.merge(bootes_df, bootesdf, on=['name'], how='left')
bootes_df = pd.merge(bootes_df, bootesbestdf, on=['name'], how='left')

#elais field

idxdesi, idxgrg, sep2d, sep3d = search_around((elaisdf['ra'], elaisdf['dec']), (elaisdesi['ra'], elaisdesi['dec']), elaisdf['z'], elaisdesi['z_phot_mean'], maxdistance)

elais_dftotal = pd.DataFrame({'name':np.array(elaisdf['name'])[idxgrg], 'ra':np.array(elaisdf['ra'])[idxgrg], 'dec':np.array(elaisdf['dec'])[idxgrg], 'z':np.array(elaisdf['z'])[idxgrg], 'radesi':elaisdesi.iloc[idxdesi, 0],\
                        'decdesi':elaisdesi.iloc[idxdesi, 1], 'sep3d': sep3d})

table = elais_dftotal['name'].value_counts()
elais_df = pd.DataFrame({'name':table.index, 'ngalaxies':table.values})
elais_df = pd.merge(elais_df, elaisdf, on=['name'], how='left')
elais_df = pd.merge(elais_df, elaisbestdf, on=['name'], how='left')


# Lockman field

idxdesi, idxgrg, sep2d, sep3d = search_around((lhdf['ra'], lhdf['dec']), (lhdesi['ra'], lhdesi['dec']), lhdf['z'], lhdesi['z_phot_mean'], maxdistance)

lh_dftotal = pd.DataFrame({'name':np.array(lhdf['name'])[idxgrg], 'ra':np.array(lhdf['ra'])[idxgrg], 'dec':np.array(lhdf['dec'])[idxgrg], 'z':np.array(lhdf['z'])[idxgrg], 'radesi':lhdesi.iloc[idxdesi, 0],\
                        'decdesi':lhdesi.iloc[idxdesi, 1], 'sep3d': sep3d})

table = lh_dftotal['name'].value_counts()
lh_df = pd.DataFrame({'name':table.index, 'ngalaxies':table.values})
lh_df = pd.merge(lh_df, lhdf, on=['name'], how='left')
lh_df = pd.merge(lh_df, lhbestdf, on=['name'], how='left')



# merge the 3 deep fields analyses 
df = pd.concat([bootes_df, elais_df, lh_df])
df = df[df.SFR_cons > -10]


df_herg = df[df['AGN_final']==1]
df_lerg = df[df['AGN_final']==0]


fig1, ax = plt.subplots(figsize = (6, 6))
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
# plt.scatter(df_herg.ngalaxies, df_herg.SFR_cons, marker='.', s = 100, edgecolors='black', facecolors='blue', label='HERG')
# plt.scatter(df_lerg.ngalaxies, df_lerg.SFR_cons, marker='.', s = 100, edgecolors='black', facecolors='orange', label='LERG')
plt.scatter(df.ngalaxies, df.power, marker='.', s = 100, edgecolors='black', facecolors='blue')

# plt.ylabel(r'$Log(SFR [M_{\odot} \ yr^{-1}]$)', fontsize=15) #20 with fraction
plt.ylabel(r'$ P_{150} [W \ Hz^{-1}] $', fontsize=15)
plt.xlabel(r'$N_{gal}$', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xscale('log')
# plt.xscale('log')
plt.yscale('log')
# ax.get_xaxis().set_major_formatter(ScalarFormatter())
# ax.get_yaxis().set_major_formatter(ScalarFormatter())
# ax.set_yticks([2e-2, 3e-2, 4e-2, 5e-2])
# ax.set_xticks([1,3,5,8,10])
plt.legend(fontsize=15)
plt.show()
exit()