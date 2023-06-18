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

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
file = main + 'BLDF/File/bootes_GRGs_catalogue.csv'
bootesbest = main + '/BLDF/File/Crossmatch/full_Bootes_Best_crossmatch.csv'

bootesdf = pd.read_csv(file, delimiter=None, usecols=('name', 'ra', 'dec', 'z', 'zstd', 'lls'))
bootesbestdf = pd.read_csv(bootesbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
bootesdf = pd.merge(bootesdf, bootesbestdf, on=['name'], how='left')


# Elais field
elaisfile = main + '/ELDF/File/en1_GRGs_catalogue.csv'
elaisbest = main + '/ELDF/File/Crossmatch/full_en1_Best_crossmatch.csv'


elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'ra', 'dec', 'z', 'zstd', 'lls'))
elaisbestdf = pd.read_csv(elaisbest, delimiter=',', usecols=('name', 'AGN_final', 'Mass_cons','SFR_cons'))
elaisdf = pd.merge(elaisdf, elaisbestdf, on=['name'], how='left')

# Lockman field 

lhfile = main + '/LHLDF/File/lh_GRGs_catalogue.csv'
lhbest = main + '/LHLDF/File/Crossmatch/full_lh_Best_crossmatch.csv'


lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'ra', 'dec', 'z', 'zstd', 'lls'))
lhbestdf = pd.read_csv(lhbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
lhdf = pd.merge(lhdf, lhbestdf, on=['name'], how='left')

df = pd.concat([bootesdf, elaisdf, lhdf])
df = df[df.z < 2]
df.SFR_cons = (df.SFR_cons)

madausfr = np.log10(np.min(df.z)*(1+df.z)**2.7)

binedges = np.linspace(-2,3, num=20)
n, _ = np.histogram(df.SFR_cons, binedges)

madaun, _ = np.histogram(madausfr, binedges)
madaun = madaun * (len(df.z) / np.sum(madaun))


plt.bar(binedges[:-1], n, width=np.diff(binedges), align='edge', alpha=0.5, label='Observed')
plt.bar(binedges[:-1], madaun, width=np.diff(binedges), align='edge', alpha=0.5, label='Expected')
plt.xlabel('SFR')
plt.ylabel('Number of Galaxies')
plt.title('Observed vs. Expected SFR Distribution')
plt.legend()
plt.show()
