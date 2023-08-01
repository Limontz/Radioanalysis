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
# from cutout2D import cutout
cosmo =  FlatLambdaCDM(H0=70, Om0=0.3)
z=0.37


##########################################################################################
from astropy import wcs

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
file = main + '/BLDF/File/Lolls_bootes_rgs.csv'

df = pd.read_csv(file)


for i in range(len(df.name)):

    spidx_file = main + '/BLDF/File/spidxmap/output/lba-hba/' + str(df.name[i]) + '_spidx.fits'

    fig1 = plt.figure(figsize = (8, 6))
    fig1.subplots_adjust(top=0.950,bottom=0.13,left=0.13,right=0.95)

    f = aplpy.FITSFigure(spidx_file, figure=fig1)

    ax = plt.gca()
    # f.axis_labels.hide()
    # f.remove_grid()
    # f.tick_labels.hide()
    # f.ticks.hide()

    f.show_colorscale(cmap='rainbow' , stretch='linear')#, vmax=max_val)#, smooth=1) #vmid=1e-8,


    # f.show_contour(lba_new_fitsfile, levels=levs, smooth=1, colors='white' , linewidths=0.4)

    ########## GRIGLIA
    #f.add_grid()

    ########## COLOR BAR
    f.add_colorbar()#ticks=[min_val, 0.0005,0.001,0.005,0.01,
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

    f.colorbar.set_axis_label_text(r'$\alpha$')
    f.colorbar.set_axis_label_font(size = 15)
    f.colorbar.set_font(size = 15)
    f.colorbar.set_location('right')
    # f.colorbar.ax.tick_params(labelsize= 15)
    f.axis_labels.set_font(size = 15)
    f.tick_labels.set_font(size= 15)
    f.savefig(main + '/BLDF/File/spidxmap/output/lba-hba/'+ str(df.name[i]) + '_spidx.png')

    # plt.show()
    plt.clf()

