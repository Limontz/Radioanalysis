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
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from astropy import units as u
import seaborn as sns
import matplotlib.image as mpimg
from PIL import Image
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter, NullFormatter
import pandas as pd
Mpc_to_m = 3.08e22
Jy_to_W = 1.e-26

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'


elaisfile = main + 'ELDF/File/en1_GRGs_catalogue.csv'
elaisfluxfile = main + 'ELDF/File/en1_grg_flux.csv'
elaiscorefile = main + 'ELDF/File/en1_grg_core_flux.csv'

bootesfile = main + 'BLDF/File/bootes_GRGs_catalogue.csv'
bootesfluxfile = main + 'BLDF/File/bootes_grg_flux.csv'
bootescorefile = main + 'BLDF/File/bootes_grg_core_flux.csv'

lhfile = main + 'LHLDF/File/LH_GRGs_catalogue.csv'
lhfluxfile = main + 'LHLDF/File/LH_grg_flux.csv'
lhcorefile = main + 'LHLDF/File/LH_grg_core_flux.csv'

elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
elaisfluxdf = pd.read_csv(elaisfluxfile, delimiter=',')
elaiscoredf = pd.read_csv(elaiscorefile, delimiter=',')

elaisdf = pd.merge(elaisdf,elaiscoredf,on='name', how='left')
elaisdf = pd.merge(elaisdf,elaisfluxdf,on='name', how='left')

bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
bootesfluxdf = pd.read_csv(bootesfluxfile, delimiter=',')
bootescoredf = pd.read_csv(bootescorefile, delimiter=',')

bootesdf = pd.merge(bootesdf,bootescoredf,on='name', how='left')
bootesdf = pd.merge(bootesdf,bootesfluxdf,on='name', how='left')

lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))
lhfluxdf = pd.read_csv(lhfluxfile, delimiter=',')
lhcoredf = pd.read_csv(lhcorefile, delimiter=',')

lhdf = pd.merge(lhdf,lhcoredf,on='name', how='left')
lhdf = pd.merge(lhdf,lhfluxdf,on='name', how='left')


df = pd.concat([elaisdf, bootesdf, lhdf])

df['corefrac'] = df.coreflux/df.flux
df['corefracerr'] = df.corefrac * ( (df.sigma_flux/df.flux)**2 + (df.fluxerr/df.coreflux)**2 )

fig = plt.figure(figsize=(6, 6))  # Adjust the figure size as needed
fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

# Normalize histograms by dividing counts by the total number of elements

bins = np.histogram((0, np.max(df.corefrac)), bins=15)[1]
sns.histplot(df.corefrac, bins=bins, kde=True, alpha=0.5, color='blue',stat='probability')

# Add labels and title
plt.xlabel(r'$f$', fontsize=15)
plt.ylabel('N', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend()
plt.show()
