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

def search_around(coord1, coord2, z1, z2, seplimit):
    c1 = SkyCoord(ra=coord1[0]*u.deg, dec=coord1[1]*u.deg, distance=cosmo.comoving_distance(z1), frame='icrs')
    # c2 = SkyCoord(ra=datadesi['ra'].head(len)*u.deg, dec=datadesi['dec'].head(len)*u.deg, distance=cosmo.comoving_distance(datadesi['z_phot_mean'].head(len)), frame='icrs')
    c2 = SkyCoord(ra=coord2[0]*u.deg, dec=coord2[1]*u.deg, distance=cosmo.comoving_distance(z2), frame='icrs')
    idx2, idx1, sep2d, sep3d = SkyCoord.search_around_3d(c1, c2, seplimit*u.Mpc)
    return(idx2, idx1, sep2d, sep3d)

#******************* file reading and data management ****************************************

# Bootes field
main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
file = main + 'BLDF/File/bootes_RGs_catalogue.csv'
desi_file = main + '/File/DESI_EDR_3df.csv'

bootesdf = pd.read_csv(file, delimiter=None, usecols=('name', 'ra', 'dec', 'z', 'ztype', 'zstd', 'lls'))
bootesdesi = pd.read_csv(desi_file, delimiter=None, usecols=('ra', 'dec', 'z_sp'))

# bootesdf = bootesdf[bootesdf.ztype == 's']

bootesdesi = bootesdesi[bootesdesi['z_sp']>0]


# dfgrg_copy = dfgrg
# datadesi_copy = datadesi

# Elais field
elaisfile = main + '/ELDF/File/en1_RGs_catalogue.csv'


elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'ra', 'dec', 'z', 'ztype', 'zstd', 'lls'))
elaisdesi = pd.read_csv(desi_file, delimiter=',', usecols=('ra', 'dec', 'z_sp'))

# elaisdf = elaisdf[elaisdf.ztype == 's']
elaisdesi = elaisdesi[elaisdesi['z_sp']>0]



depthg = 24.0
depthr = 23.4
depthz = 22.5

zmax = 0.7 #np.max(dfgrg['z'])  #
elaisdf = elaisdf[elaisdf['z'] <= zmax]
bootesdf = bootesdf[bootesdf['z'] <= zmax]


depthg = 24.0
depthr = 23.4
depthz = 22.5

# mdesi_bootes = depthr - 5 * np.log10(cosmo.luminosity_distance(zmax).to(u.pc)/u.pc) + 5 * np.log10(cosmo.luminosity_distance(bootesdesi['z_phot_mean']).to(u.pc)/u.pc) #+Kcorrection?
# mdesi_elais = depthr - 5 * np.log10(cosmo.luminosity_distance(zmax).to(u.pc)/u.pc) + 5 * np.log10(cosmo.luminosity_distance(elaisdesi['z_phot_mean']).to(u.pc)/u.pc) #+Kcorrection?

# elaisdesi = elaisdesi[elaisdesi['dered_mag_r'] <= mdesi_elais]
# bootesdesi = bootesdesi[bootesdesi['dered_mag_r'] <= mdesi_bootes]

#*********************************************************************************************************************************************************
#******************* number of galaxies around the grgs **************************************************************************************************
#*********************************************************************************************************************************************************

maxdistance = 10 #Mpc unit

# Bootes field
idxdesi, idxgrg, sep2d, sep3d = search_around((bootesdf['ra'], bootesdf['dec']), (bootesdesi['ra'], bootesdesi['dec']), bootesdf['z'], bootesdesi['z_sp'], maxdistance)

bootes_dftotal = pd.DataFrame({'name':np.array(bootesdf['name'])[idxgrg], 'ra':np.array(bootesdf['ra'])[idxgrg], 'dec':np.array(bootesdf['dec'])[idxgrg], 'z':np.array(bootesdf['z'])[idxgrg], 'radesi':bootesdesi.iloc[idxdesi, 0],\
                        'decdesi':bootesdesi.iloc[idxdesi, 1], 'sep3d': sep3d})

table = bootes_dftotal['name'].value_counts()
bootes_df = pd.DataFrame({'name':table.index, 'ngalaxies':table.values})
bootes_df = pd.merge(bootes_df, bootesdf, on=['name'], how='left')

#elais field

idxdesi, idxgrg, sep2d, sep3d = search_around((elaisdf['ra'], elaisdf['dec']), (elaisdesi['ra'], elaisdesi['dec']), elaisdf['z'], elaisdesi['z_sp'], maxdistance)

elais_dftotal = pd.DataFrame({'name':np.array(elaisdf['name'])[idxgrg], 'ra':np.array(elaisdf['ra'])[idxgrg], 'dec':np.array(elaisdf['dec'])[idxgrg], 'z':np.array(elaisdf['z'])[idxgrg], 'radesi':elaisdesi.iloc[idxdesi, 0],\
                        'decdesi':elaisdesi.iloc[idxdesi, 1], 'sep3d': sep3d})

table = elais_dftotal['name'].value_counts()
elais_df = pd.DataFrame({'name':table.index, 'ngalaxies':table.values})
elais_df = pd.merge(elais_df, elaisdf, on=['name'], how='left')


# merge the 3 deep fields analyses 
df_3df = pd.DataFrame({'name': np.concatenate((bootes_df.name, elais_df.name), axis=None),\
                       'ngalaxies': np.concatenate( (bootes_df.ngalaxies,elais_df.ngalaxies,), axis=None),\
                       'lls': np.concatenate( (bootes_df.lls,elais_df.lls), axis=None)})


df_3df.to_csv(main + '/File/photo_spec_total_ngalaxies.csv')

exit()
binwidth = 10
xymax = max(np.max(np.abs(df_3df['ngalaxies']/(np.pi*1.e2))), np.max(np.abs(df_3df['ngalaxies']/(np.pi*1.e2))))
xymin = min(np.min(np.abs(df_3df['ngalaxies']/(np.pi*1.e2))), np.min(np.abs(df_3df['ngalaxies']/(np.pi*1.e2))))


lim = (int(np.log10(xymax)/binwidth) + 1) * binwidth

bins = np.arange(-lim, lim + binwidth, binwidth)
bins = 10 ** np.linspace(np.log10(1), np.log10(xymax), 20)


fig1 = plt.figure(figsize = (6, 6))
fig1, ax = plt.subplots(figsize = (6, 6))
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
#
rg = df_3df['ngalaxies'][df_3df['lls'] < 0.5]
grg = df_3df['ngalaxies'][df_3df['lls'] >= 0.5]

print(np.mean(rg), np.mean(grg), np.max(rg), np.max(grg) )
print(ks_2samp(rg, grg))

sns.histplot(np.log10(rg/(np.pi*1.e2)), bins=15, edgecolor="white", color = 'dimgray', zorder=2, label = r'$ LLS < 0.5 Mpc$', binrange=(np.log10(xymin), np.log10(xymax)))
sns.histplot(np.log10(grg/(np.pi*1.e2)), bins=15, edgecolor="white", color = 'darkblue', zorder=3, label = r'$ LLS \geq 0.5 Mpc $', alpha = 0.5, binrange=(np.log10(xymin), np.log10(xymax)))
plt.grid(color = 'white', zorder=0, linewidth = 1.5)
plt.xlabel(r'$Log(\Sigma_{gal}  \ \rm  [Mpc^{-2}])$',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'$N(Log(\Sigma_{gal}))$', fontsize = 15) #20 with fraction
plt.yticks(fontsize=15)
plt.xticks([-1.5, -1.0, -0.5, 0.0, 0.5 ], fontsize=15)
ax.patch.set_facecolor('#ababab')
ax.patch.set_alpha(0.3)
plt.legend(fontsize=15)
plt.yscale('log')
plt.show()

exit()

#******************************************************************************
#************************************* ngal distribution **********************
#******************************************************************************
rmax = maxdistance #u.Mpc
bins = 25
bootes_dftotal = bootes_dftotal[bootes_dftotal['sep3d'] <= rmax*u.Mpc]
elais_dftotal = elais_dftotal[elais_dftotal['sep3d'] <= rmax*u.Mpc]

volume = maxdistance**3/ (bins)

r = [((i)*volume)**(1/3) for i in range(bins)]
test = [r[i+1]**3 - r[i]**3 for i in range(len(r)-1)]

# Bootes field
ngal_rg = np.zeros((len(bootes_df[bootes_df['lls']<0.7]),bins-1))
ngal_grg = np.zeros((len(bootes_df[bootes_df['lls']>=0.7]),bins-1))

m = 0
n = 0

for name in bootes_df['name']:

    if (bootes_df.loc[bootes_df['name'] == name, 'lls'].iloc[0] < 0.7):
       # if (df.loc[df['name'] == name, 'sep3d'].iloc[0] < 0.7):
       ngal_rg[m], binedges = np.histogram(bootes_dftotal['sep3d'][bootes_dftotal['name'] == name], bins=r, range=(0,rmax))
       m += 1
    else:
       ngal_grg[n], binedges = np.histogram(bootes_dftotal['sep3d'][bootes_dftotal['name'] == name], bins=r, range=(0,rmax))
       n += 1
    radius = [(t + s)/2 for s, t in zip(binedges, binedges[1:])]


#elais field

m = 0
n = 0

elais_ngal_rg = np.zeros((len(elais_df[elais_df['lls']<0.7]),bins-1))
elais_ngal_grg = np.zeros((len(elais_df[elais_df['lls']>=0.7]),bins-1))

for name in elais_df['name']:

    if (elais_df.loc[elais_df['name'] == name, 'lls'].iloc[0] < 0.7):
       elais_ngal_rg[m], binedges = np.histogram(elais_dftotal['sep3d'][elais_dftotal['name'] == name], bins=r, range=(0,rmax))
       m += 1
    else:
       elais_ngal_grg[n], binedges = np.histogram(elais_dftotal['sep3d'][elais_dftotal['name'] == name], bins=r, range=(0,rmax))
       n += 1
    radius = [(t + s)/2 for s, t in zip(binedges, binedges[1:])]

#Lockman field

m = 0
n = 0

lockman_ngal_rg = np.zeros((len(lh_df[lh_df['lls']<0.7]),bins-1))
lockman_ngal_grg = np.zeros((len(lh_df[lh_df['lls']>=0.7]),bins-1))

for name in lh_df['name']:

    if (lh_df.loc[lh_df['name'] == name, 'lls'].iloc[0] < 0.7):
       lockman_ngal_rg[m], binedges = np.histogram(lh_dftotal['sep3d'][lh_dftotal['name'] == name], bins=r, range=(0,rmax))
       m += 1
    else:
       lockman_ngal_grg[n], binedges = np.histogram(lh_dftotal['sep3d'][lh_dftotal['name'] == name], bins=r, range=(0,rmax))
       n += 1
    radius = [(t + s)/2 for s, t in zip(binedges, binedges[1:])]


ngal_grg = np.concatenate((ngal_grg, np.concatenate( (elais_ngal_grg, lockman_ngal_grg), axis=0) ), axis=0)
ngal_rg = np.concatenate((ngal_rg, np.concatenate( (elais_ngal_rg, lockman_ngal_rg), axis=0) ), axis=0)

area = np.pi*((radius[1]-radius[0]))**2
mean_rg = np.mean(ngal_rg/1.e2, axis = 0)
mean_grg = np.mean(ngal_grg/1.e2, axis = 0)
sigma_rg = np.std(ngal_rg/1.e2, axis = 0)
sigma_grg = np.std(ngal_grg/1.e2, axis = 0)

fig1 = plt.figure(figsize = (6, 6))
fig1, ax = plt.subplots(figsize = (6, 6))
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

# plt.bar(radius, mean_rg*1.e6, width=binedges[1]-binedges[0], align='edge', color='grey', edgecolor='white', linewidth=1, zorder=2, label='RGs')
# plt.bar(radius, mean_grg*1.e6, width=binedges[1]-binedges[0], align='edge', color='blue', edgecolor='white', alpha = 0.5, linewidth=1, zorder=3, label='GRGs')
# plt.grid(color = 'white', zorder=0, linewidth = 1.5)
# plt.xlabel(r'$D \ [Mpc]$',\
#            fontsize = 15) #20 with fraction
# plt.ylabel(r'$\Sigma_{gal} \ [kpc^{-2}] \cdot 10^6$', fontsize = 15) #20 with fraction
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)
# ax.patch.set_facecolor('#ababab')
# ax.patch.set_alpha(0.3)
# plt.xlim(0, 1300)
# plt.xscale('log')
# ax.get_xaxis().set_major_formatter(ScalarFormatter())
# radius[0] = 0
plt.plot(radius, mean_rg*1.e2, linestyle = '-', marker = '.', markersize = 5.0, color= 'grey', label = 'RGs', linewidth = 2.0)
plt.plot(radius, mean_grg*1.e2, linestyle = '-', marker = '.', markersize = 5.0, color= 'blue', label = 'GRGs', linewidth = 2.0)
plt.xlabel(r'$D$ [Mpc]',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'$N_{gal}$', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# plt.xscale('log')
# plt.yscale('log')
ax.get_xaxis().set_major_formatter(ScalarFormatter())
ax.get_yaxis().set_major_formatter(ScalarFormatter())
# ax.set_yticks([2e-2, 3e-2, 4e-2, 5e-2])
ax.set_xticks([1,3,5,8,10])
plt.legend(fontsize=15)
plt.show()
exit()
#***************************************************************************************************************
#******************************* orientation of jets vs galaxies distribution **********************************
#***************************************************************************************************************

def orientation(coord1, coord2): #position angle [0,180)

    c1 = SkyCoord(ra=coord1[0]*u.deg, dec=coord1[1]*u.deg, frame='icrs')
    c2 = SkyCoord(ra=coord2[0]*u.deg, dec=coord2[1]*u.deg, frame='icrs')

    pa = c1.position_angle(c2).to(u.deg)
    pa = np.where(pa > 180*u.deg, 360*u.deg - pa, pa)
    return (pa)


pa = orientation((dftotal['ra'],dftotal['dec']), (dftotal['radesi'], dftotal['decdesi']))
dftotal['padesi'] = pa

dftotal = pd.merge(dftotal, dfpa, on=['name'], how='left')
delta = (dftotal['padesi'] - dftotal['pa'])
delta = np.where(((delta >= -90*u.deg)&(delta<90*u.deg)), np.abs(delta),\
        np.where(((delta >= -180*u.deg)&(delta<-90*u.deg)), delta+180,\
        np.where(((delta >90*u.deg)&(delta<=180*u.deg)), 180-delta, delta)))


from scipy.stats import uniform
from scipy.stats import ks_2samp
from scipy.stats import kstest
from scipy import stats
# print("--- %s seconds ---" % (time.time() - start_time))
# unif = uniform.rvs(size=8000, loc=0, scale=90)
delta = delta[~np.isnan(delta)]
# for i in delta:
#     print(i)

p = []
for i in range(1000):
    unif=stats.uniform(0,90).rvs(len(delta))
    # print(np.max(unif))

# unif2=np.random.uniform(low=0, high=90, size=10000)
    p.append(ks_2samp(delta, unif)[-1])
# exit()


# fig1, ax = plt.subplots(figsize = (6, 6))
# fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
#
# # sns.histplot(unif, binwidth=5, edgecolor="white", color = 'grey', zorder=2, label = 'RGs', stat='probability')
# sns.histplot(delta, binwidth=5, edgecolor="white", color = 'darkblue', zorder=3, label = 'RGs', stat='probability')
# plt.grid(color = 'white', zorder=0, linewidth = 1.5)
# plt.xlabel(r'$\Delta \theta$', fontsize = 15) #20 with fraction
# plt.ylabel('Density',fontsize = 15) #20 with fraction
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)
# ax.patch.set_facecolor('#ababab')
# ax.patch.set_alpha(0.3)
#
# # ax2 = ax.twinx()
# # sns.histplot(delta, ax=ax2, binwidth=10, edgecolor="white", color = 'darkblue', zorder=3, label = 'RGs')
# # ax2.set_ylabel('N', fontsize=15)
# # ax2.tick_params(axis='y', labelsize=15)
# # ax.set_xscale('symlog')
# # ax.get_xaxis().set_major_formatter(ScalarFormatter())
# # plt.legend(fontsize=15)
# plt.show()
# exit()

print(df)
for name in dftotal['name']:
    # if (df.loc[df['name'] == name, 'lls'].iloc[0] >= 0.7):
       # print(dftotal['padesi'][dftotal['name'] == name], dftotal['radesi'][dftotal['name'] == name], dftotal['decdesi'][dftotal['name'] == name])

       ax.scatter(217.2997,33.4437\
                  , marker = [2, 2, 0],\
                  s = 1000 * df['lls'][df['name'] == '1429+3326'], color = 'red')
       ax.scatter(dftotal['radesi'][dftotal['name'] == '1429+3326'], dftotal['decdesi'][dftotal['name'] == '1429+3326'],\
                  marker = '.', s = 20.0, color = 'black')
       plt.xlabel('RA (J2000)', fontsize=15)
       plt.ylabel('Dec (J2000)', fontsize=15)
       plt.xticks(fontsize=15)
       plt.yticks(fontsize=15)
       ax.invert_xaxis()
       plt.show()
       plt.clf()
