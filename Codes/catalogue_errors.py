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

file = main + 'BLDF/File/bootes_RGs_catalogue.csv'

df = pd.read_csv(file, delimiter=',')
# print(np.mean(df.z))
# exit()
df = df[df.z==1.5]
print(len(df.name))
exit()
# look for doubles
print('duplicate RG index: ', np.where(df.duplicated(subset=['name'], keep='first'))[0])
exit()
# name and coordinates match 
new_list = [x[0:4] + x[9:14] for x in df.coords ]
# print('RG index with coords-name mismatch: ',np.where(~(df.name == np.array(new_list)))[0])

#type of z
# print('type of redshift: ', [*set(df.ztype)])


#type of survey: should be either SDSS, DESI or CWISE
# print('survey: ', [*set(df.survey)])

# survey-Hostname mismatch
new_list = np.array([x for x in df.Hostname[df.survey=='SDSS']])
index = np.where(~(np.array([x[0:2] == 'J1' for x in new_list])))
# print( 'SDSS Hostname: ', new_list[index])


new_list = np.array([x for x in df.Hostname[df.survey=='CWISE']])
index = np.where(~(np.array([x[0:2] == 'J1' for x in new_list])))
# print( 'CWISE Hostname: ', new_list[index])

new_list = np.array([x for x in df.Hostname[df.survey=='DESI']])
index = np.where(~(np.array([x[0:3] == 'J16' for x in new_list])))
# print( 'DESI Hostname: ', new_list[index])


# galaxy type
# print('galaxy type: ', [*set(df.type)])


# survey-magnitude mismatch
df.loc[(df["survey"] == "SDSS") & (df["mag"].str[-1] == "r"), "mag"] = df["mag"].str[:-1] + "r'"

new_list = np.array([x for x in df.mag[df.survey == 'SDSS']])
index = np.where(~(np.array([x[5:7] == "r'" for x in new_list])))
print( 'SDSS magnitude: ', new_list[index])

new_list = np.array([x for x in df.mag[df.survey == 'CWISE']])
index = np.where(~(np.array([x[5:7] == "W1" for x in new_list])))
print( 'CWISE magnitude: ', new_list[index])


new_list = np.array([x for x in df.mag[df.survey == 'DESI']])
index = np.where(~(np.array([x[5:7] == "r" for x in new_list])))
print( 'DESI magnitude: ', new_list[index])

# FR type
print('galaxy type: ', [*set(df.FRtype)])
exit()