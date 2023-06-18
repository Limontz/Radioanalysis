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

bootesfile = main + 'BLDF/File/bootes_GRGs_catalogue.csv'
bootesfluxfile = main + 'BLDF/File/bootes_grg_flux.csv'

lhfile = main + 'LHLDF/File/LH_GRGs_catalogue.csv'
lhfluxfile = main + 'LHLDF/File/LH_grg_flux.csv'

elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
elaisfluxdf = pd.read_csv(elaisfluxfile, delimiter=',')
elaisdf = pd.merge(elaisdf,elaisfluxdf,on='name', how='left')

bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
bootesfluxdf = pd.read_csv(bootesfluxfile, delimiter=',')
bootesdf = pd.merge(bootesdf,bootesfluxdf,on='name', how='left')

lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))
lhfluxdf = pd.read_csv(lhfluxfile, delimiter=',')
lhdf = pd.merge(lhdf,lhfluxdf,on='name', how='left')

df = pd.concat([elaisdf, bootesdf, lhdf])

alpha = 0.7
df['power'] = 4. * np.pi * (1 + df.z)**(alpha - 1) * (df.flux/1e3) *  ((cosmo.luminosity_distance(df.z) * Mpc_to_m / u.Mpc)**2) * Jy_to_W   #units W


Dabhade_2020 = main + 'BLDF/File/Dabhade(2020)sample.grg'
Dabhade_2017 = main + 'BLDF/File/Dabhade(2017)sample.grg'
Kuzmicz = main + 'BLDF/File/Kuzmicz(2021)sample.grg'
Heinz_file = main + 'BLDF/File/178GRGs.dat'
sagan_file = main + 'BLDF/File/SAGANI_grg.csv'
bhukta_file = main + 'BLDF/File/Bhukta_tgss_grg.csv'

sagan = pd.read_csv(sagan_file, delimiter=',', usecols=(10,18))
# print(sagan)
# exit()
# bhukta = pd.read_csv(bhukta_file, delimiter=',', usecols=(9,14))

vla_freq = 14200
lofar_freq = 150
spdix = 0.7

with open(Heinz_file, 'r') as data:
        H2021_power = []
        H2021_LLS = []
        for line in data:
            if not line.startswith("#"):
               p = line.split()
               H2021_power.append(10**(float(p[14])))
               # z.append(float(p[5]))
               H2021_LLS.append(float(p[8]))


H2021_power = np.array(H2021_power)
H2021_power = H2021_power * (888/150)**0.7
H2021_LLS = np.array(H2021_LLS)


with open(Dabhade_2020, 'r') as data:
        D2020_power = []
        D2020_error_power = []
        D2020_LLS = []
        for line in data:
            if not line.startswith("#"):
               p = line.split()
               D2020_power.append(float(p[16]) * 1e26)
               D2020_error_power.append(float(p[18]) *1e26)
               D2020_LLS.append(float(p[12]))

D2020_power = np.array(D2020_power)
D2020_error_power = np.array(D2020_error_power)
D2020_LLS = np.array(D2020_LLS)


with open(Dabhade_2017, 'r') as data:
        D2017_power = []
        D2017_error_power = []
        D2017_LLS = []
        for line in data:
            if not line.startswith("#"):
               p = line.split()
               D2017_power.append(float(p[11]) * 1e25)
               D2017_error_power.append(float(p[13]) * 1e25)
               D2017_LLS.append(float(p[7]))

D2017_power = np.array(D2017_power)
D2017_power = D2017_power * (1400/150)**0.7
D2017_error_power = np.array(D2017_error_power)
D2017_LLS = np.array(D2017_LLS)


with open(Kuzmicz, 'r') as data:
        Kuzmicz_power = []
        Kuzmicz_LLS = []
        for line in data:
            if not line.startswith("#"):
               p = line.split()
               Kuzmicz_power.append(10**(float(p[12])))
               Kuzmicz_LLS.append(float(p[9]))


Kuzmicz_power = np.array(Kuzmicz_power)
Kuzmicz_power = Kuzmicz_power * (1400/150)**0.7
Kuzmicz_LLS = np.array(Kuzmicz_LLS)


cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
alpha = 0.7

# power = 4. * np.pi * (1 + z)**(alpha - 1) * (cosmo.luminosity_distance(z) * Mpc_to_m / u.Mpc)**2\
        # * (flux / 1.e3) * Jy_to_W #units W


fig1, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, sharey=True, figsize=(7,7))
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0,wspace=0)
ax1.plot(sagan.iloc[:, 0], sagan.iloc[:, 1]*1.e25,linestyle = '', marker = "1", markersize = 8.0, label = 'Dabhade(2020b)')
# plt.plot(sagan.iloc[:, 0], bhukta.iloc[:, 1]*1.e25,linestyle = '', marker = "P", markersize = 4.0, label = 'Bhukta(2021)')
ax1.plot(D2017_LLS, D2017_power, linestyle = '', marker = '*', markersize = 8.0, label = 'Dabhade(2017)' )
ax1.plot(D2020_LLS, D2020_power, linestyle = '', marker = 's', markersize = 3.0, label = 'Dabhade(2020a)')
ax1.plot(Kuzmicz_LLS, Kuzmicz_power, linestyle = '', marker = 'D', markersize = 3.0, label = 'Kuzmicz(2021)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_yticklabels(labels='', fontsize=15)
ax1.legend(loc = 'lower right')


ax2.plot(H2021_LLS, H2021_power, linestyle = '', marker = 'X', markersize = 4.0, label = 'Andernach(2021)')
ax2.plot(df.lls, df.power, linestyle = '', marker = '.', markersize = 10.0, label = r'DF-GRGs', color = 'black')
ax2.plot(df.lls, df.power, linestyle = '', marker = '.', markersize = 10.0, color = 'black')

ax2.set_ylabel(r'$P_{150 MHz} \rm [W \ Hz^{-1}]$', fontsize=15) #20 with fraction
ax2.yaxis.set_label_coords(-0.1,1)
plt.xlabel(r'$LLS$ [Mpc]', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
plt.xscale('log')
ax1.set_xticks([0.7, 1, 2, 3, 4])
ax1.get_xaxis().set_major_formatter(ScalarFormatter())
plt.legend(loc = 'lower right')
plt.show()
