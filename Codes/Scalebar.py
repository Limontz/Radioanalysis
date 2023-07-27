def addScalebar(ax, cdelt, z, kpc,  fontsize, color='white'):
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from astropy.cosmology import FlatLambdaCDM
    import numpy as np
    import astropy.units as u
    # cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    # size = 1/wcs
    # print("-- Redshift: %f" % z)
    # degperpixel = np.abs(wcs.pixel_to_world(0,0).dec - wcs.pixel_to_world(0,1).dec) # delta deg for 1 pixel
    # degperkpc = cosmo.arcsec_per_kpc_proper(z).value/3600.
    # pixelperkpc = degperkpc/degperpixel
    fontprops = fm.FontProperties(size=fontsize)
    size = (1/60)/cdelt
    # print(size)
    scalebar = AnchoredSizeBar(ax.transData, size, r'$1^\prime$','lower right', fontproperties=fontprops, pad=0.5, color=color, frameon=False, sep=1, label_top=True)
    ax.add_artist(scalebar)
