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
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter, NullFormatter
from lib_linearfit import linear_fit, linear_fit_bootstrap
import pandas as pd

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
bootesfile = main + 'BLDF/File/bootes_RGs_catalogue.csv'
elaisfile = main + 'ELDF/File/en1_RGs_catalogue.csv'
lhfile = main + 'LHLDF/File/lh_RGs_catalogue.csv'

dfbootes = pd.read_csv(bootesfile, delimiter=',')
dfelais = pd.read_csv(elaisfile, delimiter=',')
dflh = pd.read_csv(lhfile, delimiter=',')

df = pd.concat([dfbootes, dfelais, dflh])

# df = df[df.ztype == 's']
# print(df)
print(np.min(df.z), np.max(df.z))


median_z = np.median(df.z)
print("median z:",median_z)
exit()
print(len(dfbootes[dfbootes.lls>=1]),len(dfelais[dfelais.lls>=1]))
x = [median_z, median_z]
y = [0,650]

binwidth = 0.20
xymax = max(np.max(np.abs(df.lls)), np.max(np.abs(df.lls)))
lim = (int(xymax/binwidth) + 1) * binwidth
bins = np.arange(0, lim + binwidth, binwidth)


fig = plt.figure(figsize = (6, 6))
fig1, ax = plt.subplots()
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
plt.plot(x, y, color = 'orange', zorder = 4, linewidth= 5.0)
sns.histplot(df.lls, bins = bins, edgecolor="white", color = 'darkblue', zorder=2)

df = df[df.ztype=='s']

# print(len(df.z))
# sns.histplot(df.z, bins = bins, edgecolor="white", color = 'darkred', zorder=3)
plt.grid(color = 'white', zorder=0, linewidth = 1.5)
plt.xlabel(r'LLS [Mpc]',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'N', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
ax.patch.set_facecolor('#ababab')
ax.patch.set_alpha(0.3)
# plt.xlim(-3.0, 4.0)
# plt.legend(fontsize=15)
# plt.xscale('log')
# plt.yscale('log')
plt.show()
