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
import math

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/File/total_ngalaxies.csv'

df = pd.read_csv(file)

binwidth = 10
xymax = max(np.max(np.abs(df['ngalaxies']/(np.pi*1.e2))), np.max(np.abs(df['ngalaxies']/(np.pi*1.e2))))
xymin = min(np.min(np.abs(df['ngalaxies']/(np.pi*1.e2))), np.min(np.abs(df['ngalaxies']/(np.pi*1.e2))))


lim = (int(np.log10(xymax)/binwidth) + 1) * binwidth

bins = np.arange(-lim, lim + binwidth, binwidth)
bins = 10 ** np.linspace(np.log10(1), np.log10(xymax), 20)


fig1 = plt.figure(figsize = (6, 6))
fig1, ax = plt.subplots(figsize = (6, 6))
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
#
rg = df['ngalaxies'][df['lls'] < 0.5]
grg = df['ngalaxies'][df['lls'] >= 0.5]

print(np.mean(rg), np.mean(grg), np.max(rg), np.max(grg) )
print(ks_2samp(rg, grg))
tick = [-1.5, -1.0, -0.5, 0.0, 0.5 ]

sns.histplot(np.log10(rg/(np.pi*1.e2)), bins=15, edgecolor="white", color = 'dimgray', zorder=2, label = r'$ LLS < 0.5 Mpc$', binrange=(np.log10(xymin), np.log10(xymax)))
sns.histplot(np.log10(grg/(np.pi*1.e2)), bins=15, edgecolor="white", color = 'darkblue', zorder=3, label = r'$ LLS \geq 0.5 Mpc $', alpha = 0.5, binrange=(np.log10(xymin), np.log10(xymax)))
plt.grid(color = 'white', zorder=0, linewidth = 1.5)
plt.xlabel(r'$Log(\Sigma_{gal}  \ \rm  [Mpc^{-2}])$',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'$N(Log(\Sigma_{gal}))$', fontsize = 15) #20 with fraction
plt.yticks(fontsize=15)
plt.xticks(tick, fontsize=15)
ax.patch.set_facecolor('#ababab')
ax.patch.set_alpha(0.3)
plt.legend(fontsize=15)
plt.yscale('log')

tick2 = np.around(10**(np.array(tick))*np.pi*1.e2, -1)
tick2[[2,3,4]] = np.around(tick2[[2,3,4]], -2)
print(tick2)

ax2 = ax.twiny()
# ax2.hist(np.log10(rg), bins=15) # Create a dummy plot
# ax2.cla()
# ax2.set_xlim(0, 4.5)
# ax2.set_xticks(fontsize=15)
ax2.set_xticklabels(np.array(tick2).astype(int), fontsize = 15)
ax2.set_xlabel('Ngal', fontsize =15)

# ax2.set_xscale('log')
plt.show()


