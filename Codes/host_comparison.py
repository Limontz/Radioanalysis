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
import pandas as pd

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'

elais_file = main + 'ELDF/file/Crossmatch/full_en1_Kondapally_crossmatch.csv'
lh_file = main + 'LHLDF/file/Crossmatch/LH_Kondapally_crossmatch.csv'

elais_duncan =  main + 'ELDF/file/Crossmatch/full_en1_Duncan_crossmatch.csv'
lh_duncan =  main + 'LHLDF/file/Crossmatch/LH_Duncan_crossmatch.csv'


elais_df = pd.read_csv(elais_file, delimiter=',')
lh_df = pd.read_csv(lh_file, delimiter=',')
elais_duncan_df = pd.read_csv(elais_duncan, delimiter=',')
lh_duncan_df = pd.read_csv(lh_duncan, delimiter=',')

elais_df= pd.merge(elais_df, elais_duncan_df, on=['NAME'], how='left')
lh_df= pd.merge(lh_df, lh_duncan_df, on=['NAME'], how='left')

elais_z = np.round(elais_df.zduncan, 3) == np.round(elais_df.Z_BEST,3)
lh_z = np.round(lh_df.zduncan,3), np.round(lh_df.Z_BEST,3)

# print( 1 - (np.count_nonzero(elais_z) + np.count_nonzero(lh_z)) / (len(elais_df.zduncan) + len(lh_df.zduncan)) )

print( 1 - (np.count_nonzero(np.round(elais_df.zduncan, 3) == np.round(elais_df.Z_BEST,3)))/len(elais_df.zduncan) )
print( 1 - (np.count_nonzero(np.round(lh_df.zduncan, 3) == np.round(lh_df.Z_BEST,3)))/len(lh_df.zduncan) )

x = [0,6]
y = [0,6]


# 2D separation

# c1 = SkyCoord(elais_df.myra*u.deg, elais_df.mydec*u.deg, frame='fk5')
# c2 = SkyCoord(elais_df.optra*u.deg, elais_df.optdec*u.deg, frame='fk5')
# sep = c1.separation(c2)
# print(len(sep))
# exit()

plt.plot(elais_df.zduncan, elais_df.Z_BEST, marker = '.', linestyle= '', color='blue', label= 'K(2021)')
# plt.plot(elais_df.las, df.LGZ_Size/60, marker = '.', linestyle= '', color='orange')
plt.plot(lh_df.zduncan, lh_df.Z_BEST, marker = '.', linestyle= '', color='blue')
# plt.plot(lh_df.las, lh_df.LGZ_Size/60, marker = '.', linestyle= '', color='orange')
plt.plot(x,y)
plt.xlabel(r'z', fontsize = 15) #20 with fraction
plt.ylabel(r'$z_K$', fontsize = 15) 
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.show()

