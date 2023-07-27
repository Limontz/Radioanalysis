from function import flatten
from astropy.wcs import WCS as pywcs
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.nddata import Cutout2D
import astropy.io.fits as fits
from astropy import units as u

def cutout(fitsfile, ra, dec, las, filename):
    header_surv, data_surv = flatten(fitsfile)
    w = pywcs(header_surv)
    size = las/header_surv['CDELT2']
    coord_source = SkyCoord(ra, dec, frame='fk5')
    data_cut = Cutout2D(data_surv, position=coord_source, size=[size,size], wcs=w)
    # Cropped WCS
    wcs_cropped = data_cut.wcs
    # Update WCS in header
    header_surv.update(wcs_cropped.to_header())
    w = pywcs(header_surv)
    
   
    # Add comment to header
    # hdr['COMMENT'] = "= Cropped fits file ({}).".format(datetime.date.today())
    # Write cropped frame to new fits file.
    # fits.writeto('crop.fits', hdu_crop.data, hdr)
    # header_surv['CRPIX1'] = size/2
    # header_surv['CRPIX2'] = size/2
    # sky = w.pixel_to_world(header_surv['CRPIX1'], header_surv['CRPIX1'])
    # print(sky)
    # exit()
    cutout = fits.PrimaryHDU(data=data_cut.data, header=header_surv)
    file = filename #'/Users/marco/Desktop/The_Master/PhD/GRG Project/' + str(field) + '/File/' + filename + obs + '.fits'
    
    cutout.writeto(file, overwrite=True)
    return(file)
