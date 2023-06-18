import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

import numpy as  np
import pandas as pd
import matplotlib.font_manager

from uncertainties import unumpy
from scipy import optimize

import aplpy
import pyregion
from astropy.coordinates import SkyCoord
from regions import PixCoord
from regions import read_ds9
#from astropy.cosmology import Planck13 as cosmo
from astropy import wcs
from astropy.io import fits
from astropy import constants as c
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from mpl_toolkits.axes_grid1 import make_axes_locatable

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
# #********************* ASTROPY METHOD *******************************

# filename = get_pkg_data_filename(main + '/Codes/1426+3222_spidx.fits')
#
# hdu = fits.open(filename)[0]
# wcs = WCS(hdu.header)
#
# fig = plt.figure(figsize = (6, 8))
# plt.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
#
# ax = plt.subplot(projection=wcs[0,:])
# plt.imshow(hdu.data, origin='lower')
# plt.grid(color='white', ls='solid')
# plt.xlabel('Galactic Longitude')
# plt.ylabel('Galactic Latitude')
#
#
# ax = plt.subplot(projection=wcs[0,:], label='overlays')
#
# im = ax.imshow(hdu.data, origin='lower')
#
# ax.coords.grid(True, color='white', ls='solid')
# ax.coords[0].set_axislabel('Galactic Longitude')
# ax.coords[1].set_axislabel('Galactic Latitude')
#
# overlay = ax.get_coords_overlay('fk5')
# overlay.grid(color='white', ls='dotted')
# overlay[0].set_axislabel('Right Ascension (J2000)')
# overlay[1].set_axislabel('Declination (J2000)')
# #plt.xlim(216, 217)
# cbar = plt.colorbar(im, orientation="horizontal")
#
# cbar.set_label(r'spectral index', fontsize = 15)
# cbar.ax.tick_params(labelsize= 15)
# plt.show()

#******************** APL METHOD ******************************************
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

import numpy as  np
import pandas as pd
import matplotlib.font_manager

from uncertainties import unumpy
from scipy import optimize
from Scalebar import addScalebar
import aplpy
import pyregion
from astropy.coordinates import SkyCoord
from regions import PixCoord
from regions import read_ds9
#from astropy.cosmology import Planck13 as cosmo
from astropy import wcs
from astropy.io import fits
from astropy import constants as c
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from matplotlib_scalebar.scalebar import ScaleBar
from cutout2D import cutout
from function import flatten
from function import SampleImage
cosmo =  FlatLambdaCDM(H0=70, Om0=0.3)



##########################################################################################
from astropy import wcs
#Immagine radio 4D --> 2D
main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
grgfile = main + '/File/grgtable2.csv'
hba_fitsfile = main + '/File/bootes_radio_image.fits'
hba_name = 'bootes_radio_image.fits'
lba_fitsfile = main + '/File/spidx/Bootes_LoLLS.fits'
lba_name = 'Bootes_LoLLS.fits'

#name='briggs_-0.25_taper15_masked_SUB2_restbeam25_image_9_-MFS-image.fits'

dfgrg = pd.read_csv(grgfile, delimiter=",", usecols=(0,1,2,3,5,11))
print(dfgrg)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 4D
# HBA MAPS
hdu4D = fits.open(hba_fitsfile)
hdu4D.info()


data = hdu4D[0].data.squeeze() # drops the size-1 axes
data = data * 1.e3
header = hdu4D[0].header
#print(header)
mywcs = wcs.WCS(header).celestial

new_header = mywcs.to_header()

######### LOFAR
# new_header['BMAJ']= (header['BMAJ'])
# new_header['BMIN']= (header['BMIN'])
# new_header['BPA' ] = (header['BPA'])
# new_header['TELESCOP']= (header['TELESCOP'])
# new_header['ORIGIN' ] = (header['ORIGIN'])

######## VLA
#new_header['BMAJ']= float(header['HISTORY'][-18][13:])/3600
#new_header['BMIN']= float(header['HISTORY'][-17][13:])/3600
#new_header['BPA' ] = float(header['HISTORY'][-3][56:])

new_fh = fits.PrimaryHDU(data=data, header=new_header)
new_fitsfile='2D-'+hba_name
new_fh.writeto(new_fitsfile, overwrite=True)

# print('beam size [arcsec x arcsec]',header['BMAJ']*3600,header['BMIN']*3600)

#print(new_fitsfile)
print('----------------------------------------------------------------------------------------------')

lba_hdu4D = fits.open(lba_fitsfile)
lba_hdu4D.info()


lba_data = lba_hdu4D[0].data.squeeze() # drops the size-1 axes
lba_data = lba_data * 1.e3
lba_header = lba_hdu4D[0].header
#print(header)
lba_mywcs = wcs.WCS(lba_header).celestial

lba_new_header = lba_mywcs.to_header()

lba_new_fh = fits.PrimaryHDU(data=lba_data, header=lba_new_header)
lba_new_fitsfile='2D-'+lba_name
lba_new_fh.writeto(lba_new_fitsfile, overwrite=True)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 2D
hdu_list_2D = fits.open(new_fitsfile)
hdu_list_2D.info()
hdu_2D = hdu_list_2D[0]
image_data_2D = hdu_2D.data

print('----------------------------------------------------------------------------------------------')


hdu_2D.header

w = wcs.WCS(hdu_2D.header)
sky = w.pixel_to_world(160, 160)
#x, y = w.world_to_pixel(sky)
#print(x, y)

print('----------------------------------------------------------------------------------------------')

lba_hdu_list_2D = fits.open(lba_new_fitsfile)
lba_hdu_list_2D.info()
lba_hdu_2D = lba_hdu_list_2D[0]
lba_image_data_2D = lba_hdu_2D.data

print('----------------------------------------------------------------------------------------------')


lba_hdu_2D.header

lba_w = wcs.WCS(lba_hdu_2D.header)
lba_sky = lba_w.pixel_to_world(160, 160)

#################################################
#     calc noise of LoLLS image
#################################################

from function import SampleImage

map = SampleImage(lba_fitsfile)
rms = map.calc_noise()

####################################
#   set parameters for the figure
####################################
levs=[3*rms*1.e3]


########## SCALE BAR
z=0.37
kpc_amin =cosmo.kpc_proper_per_arcmin(z)
kpc_asec = kpc_amin.to(u.kpc/u.arcsec) #kpc/asec
scale_500kpc = (500*u.kpc/kpc_amin).to(u.deg)
scale_100kpc = (500*u.kpc/kpc_amin).to(u.deg)
scale_100kpc, (500*u.kpc/kpc_amin)/60

# f.add_scalebar(scale_500kpc, '500kpc', color='white') #corner='top right'
#f.scalebar.set(weight='bold',linewidth=3)
#f.scalebar.set_frame(True)

########## BEAM
# f.add_beam()
# #f.beam.set_color('lightblue')
# f.beam.set_frame(True)
#f.beam.set_hatch('+')
#f.beam.set_linestyle('dotted')

# text=str(kpc_amin.value)[:4]
# f.add_label(0.5,0.97, " 1' / "+text+" kpc",relative=True, bbox=dict(facecolor='white', edgecolor='black',  pad=5.0))

# rms_muJy=rms*1e+6
# '%.2e' % rms_muJy
# f.add_label(0.85,0.97, ' '+ '%.2e' % rms_Jy +r' Jy/beam ',relative=True, color='white')
#boxstyle='round'

########## CONTORNI
# color_list=['red', 'lime', 'lime',
#             'lime','lime','lime',
#             'lime','lime', 'lime',
#             'lime','lime', 'lime',
#             'lime']

#path2='./2D-briggs_-0.25_taper15_masked_SUB2_restbeam25_PSZ2G226.18+76.79_image_9_-MFS-image.fits'
# f.show_contour(new_fitsfile, levels=levs_neg ,smooth=1, colors=color_list )
min_val=0.05
# f.show_colorscale(cmap='gnuplot2' , stretch='linear', vmin=min_val)#, vmax=max_val)#, smooth=1) #vmid=1e-8,
# f.show_contour(lba_new_fitsfile, levels=levs, smooth=1, colors='white' , linewidths=0.25)

# fig1, ax = plt.subplots(1)
# fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

fig = plt.figure(figsize=(4, 5))
# fig, ax = plt.subplots(ncols=5, nrows=5, figsize=(6,6))
start = 0
for i in range(start,35):


    grgcutout, cdelt = cutout(hba_fitsfile, dfgrg['ra'][i]*u.deg, dfgrg['dec'][i]*u.deg, (dfgrg['las'][i] + 3.5)/60)
    lba_cutout, lba_cdelt = cutout(lba_fitsfile, dfgrg['ra'][i]*u.deg, dfgrg['dec'][i]*u.deg, (dfgrg['las'][i] + 100)/60)

    data = grgcutout.data[~np.isnan(grgcutout.data)]
    oldrms = 1000
    eps = np.nanstd(data)*1e-3
    for m in range(1000):
        rms = np.nanstd(data)
        if np.abs(oldrms-rms)/rms < eps:
           noise = rms
            #print('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))
            # logging.debug('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))

        data = data[np.abs(data)<5*rms]
        oldrms = rms
    print(noise)
    z=dfgrg['z'][i]
    # print(grgcutout.data)
    # print(dfgrg['name'][i])

    # ax = plt.subplot(2, 2, i -4 + 1)

    # sky = SkyCoord(dfgrg['ra'][i], dfgrg['dec'][i], unit='deg')
    # y, x = lba_w.world_to_pixel(sky)
    #
    # y = int(y)
    # x = int(x)
    # print(new_fitsfile)

    # z=dfgrg['z'][i]
    kpc_amin =cosmo.kpc_proper_per_arcmin(z)
    scale_500kpc = (500*u.kpc/kpc_amin).to(u.deg)
    scale_100kpc = (100*u.kpc/kpc_amin).to(u.deg)
    # scale_500kpc, (500*u.kpc/kpc_amin)/60

    f = aplpy.FITSFigure(grgcutout, figure=fig, subplot=(7,5,i-start+1))

    # f.show_colorscale(cmap='gnuplot2' , stretch='linear', vmin=noise/10, vmax=0.80*np.max(data))
    f.show_colorscale(cmap='gnuplot2' , stretch='sqrt', vmin= 0)

    # f.show_contour(lba_new_fitsfile, levels=levs, smooth=1, colors='black' , linewidths=0.8)

    ax = plt.gca()
    font = 4
    name= str(dfgrg['name'][i])
    lls= str(dfgrg['lls'][i])
    addScalebar(ax, cdelt, dfgrg['z'][i], dfgrg['lls'][i], font , color='white')
    # f.add_label(0.37,0.93, name,relative=True, size = font, color='white')
    # f.add_label(0.89,0.93, lls,relative=True, size = font, color='white')
    f.show_markers(dfgrg['ra'][i]*u.deg, dfgrg['dec'][i]*u.deg, c = '', marker= 'o', s = 20, facecolors='none', edgecolors='cyan')
    plt.subplots_adjust(wspace=0.05, hspace=0.05)

    ax = plt.gca()
    f.axis_labels.hide()
    # f.remove_grid()
    f.tick_labels.hide()
    f.ticks.hide()

    # f.axis_labels.set_font(size = 15)
    # f.tick_labels.set_font(size= 15)

########## ZOOM
    # f.recenter(dfgrg['ra'][i], dfgrg['dec'][i], width=((dfgrg['las'][i] + 1.5)/60), height=((dfgrg['las'][i] + 1.5)/60)) #219.2216  34.2807 #217.4277, 33.9486
########## REGIONS
# path_region = '/home/utente/Documents/TesiProjectHeritage/Radio/LOFAR/A1413/PROVA_ds9/beam15/Allregion_ds9_justMH_smaller.reg'
# f.show_regions(path_region)

########## GRIGLIA
#f.add_grid()

########## COLOR BAR
# f.add_colorbar()#ticks=[min_val, 0.0005,0.001,0.005,0.01,
                       #min_val*1,
                       #min_val*10,
                       #min_val+100,
                       #min_val*1000,
                       #min_val*10000,
                       #0.005,
                       #max_val],
                #log_format = True
              #)
              #axis_label_text=' Jy / beam')

# f.colorbar.set_axis_label_text('mJy/beam')
# f.colorbar.set_axis_label_font(size = 15)
# f.colorbar.set_font(size = 15)
# f.colorbar.set_location('right')
# f.colorbar.ax.tick_params(labelsize= 15)
# f.ylabel(fontisize = 15)
#f.savefig('/home/utente/Documents/TesiProjectHeritage/Radio/LOFAR/A1413/image_png/A1413_-0.25_35"_zoom.png')
# f.savefig('/home/utente/Documents/TesiProjectHeritage/Radio/LOFAR/A1413/image_png/A1413_radio_rest35_and_GRID_MH.png')
# fig.tight_layout()
plt.show()
# f.save('test_spidx.png')
