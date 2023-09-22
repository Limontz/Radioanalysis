#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os, sys, argparse, logging
import numpy as np
from lib_linearfit import linear_fit_bootstrap, linsq_spidx, twopoint_spidx_bootstrap
import lib_fits
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion
# https://github.com/astrofrog/reproject
from reproject import reproject_interp, reproject_exact
reproj = reproject_exact
logging.root.setLevel(logging.DEBUG)
import pandas as pd
from function import calculate_flux

parser = argparse.ArgumentParser(description='Make spectral index maps, e.g. spidxmap.py --region ds9.reg --noise --sigma 5 --save *fits')
# parser.add_argument('images', nargs='+', help='List of images to use for spidx')
parser.add_argument('--beam', dest='beam', nargs='+', type=float, help='3 parameters final beam to convolve all images (BMAJ (arcsec), BMIN (arcsec), BPA (deg))')
parser.add_argument('--bgreg', dest='bgreg', help='DS9 region file for background estimation.')
parser.add_argument('--region', dest='region', help='Ds9 region to restrict analysis')
# parser.add_argument('--size', dest='size', nargs=2, type=float, help='Size (horizontal and vertical) of final image in degree')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, help='RA/DEC where to center final image in deg (if not given, center on first image)')
parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before calculating spidx')
parser.add_argument('--noise', dest='noise', action='store_true', help='Calculate noise of each image')
parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate results')
parser.add_argument('--force', dest='force', action='store_true', help='Force remake intermediate results')
parser.add_argument('--sigma', dest='sigma', type=float, help='Restrict to pixels above this sigma in all images')
parser.add_argument('--fluxscaleerr', dest='fluxscaleerr', type=float, default=0.0, help='Systematic error of flux scale. One value for all images.')
parser.add_argument('--upperlimit', dest='upperlimit', type=float, help='Place upper limits below this value if not detected at highest frequency at sigma. Float, e.g. -1.0')
parser.add_argument('--lowerlimit', dest='lowerlimit', type=float, help='Place lower limits below this value if not detected at lowest frequency at sigma. Float, e.g. -0.6')
parser.add_argument('--circbeam', dest='circbeam', action='store_true', help='Force final beam to be circular (default: False, use minimum common beam area)')
# parser.add_argument('--output', dest='output', default='spidx.fits', help='Name of output mosaic (default: spidx.fits)')

args = parser.parse_args()

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/BLDF/File/'
rgfile = main + 'bootes_RGs_catalogue.csv'
file = main + 'Lolls_bootes_RGs.csv'

dfgrg = pd.read_csv(rgfile, delimiter=",", usecols=('name', 'ra', 'dec', 'las', 'z', 'lls'))
df = pd.read_csv(file)
# idx = [5,14,19,31]
# df = df.iloc[idx]
df = pd.merge(df['name'], dfgrg[['name','ra', 'dec', 'las', 'z', 'lls']], on='name')
badimage = ['1443+3606', '1432+3608', '1428+3558', '1443+3555','1431+3532', '1436+3359','1426+3339','1430+3333',
            '1433+3209','1436+3144','1435+3705']

if __name__ == '__main__':
   spidx = []
   err = []
   for i in range(42,len(df.name)):
    
    
        images = [main+'spidxmap/hba/'+str(df.name[i])+'_hba.fits', main+'spidxmap/lba/'+str(df.name[i])+'_lba.fits'] 
      #   output = main+ 'spidxmap/output/lba-hba/images/' + str(df.name[i]) + '.fits'
        radec = [df.ra[i], df.dec[i]]
        region = main + 'spidxmap/ds9regions/lba' +str(df.name[i]) + '.reg'
        size = [df.las[i]*3/60, df.las[i]*3/60]
      #   if np.any(df.name[i] == np.array(badimage)):
      #      size = [df.las[i]*2/60, df.las[i]*2/60]
      #   check input
        if len(images) < 2:
          logging.error('Requires at lest 2 images.')
          sys.exit(1)
        elif len(images) > 2 and (args.upperlimit or args.lowerlimit):
          logging.error('Upper-limit currently only supported for two frequencies')
          sys.exit()

        if args.beam is not None and len(args.beam) != 3:
          logging.error('Beam must be in the form of "BMAJ BMIN BPA" (3 floats).')
          sys.exit(1)

        if radec is not None and len(radec) != 2:
          logging.error('--radec must be in the form of "RA DEC" (2 floats).')
          sys.exit(1)


        ########################################################
        # prepare images and make catalogues if necessary
        all_images = lib_fits.AllImages([imagefile for imagefile in images])

        #####################################################
        # find+apply shift w.r.t. first image
        if args.shift:
           if all_images.suffix_exists('si-shift') and not args.force:
               logging.info('Reuse si-shift images.')
               all_images = lib_fits.AllImages([name.replace('.fits', '-si-shift.fits') for name in all_images.filenames])
           else:
               all_images.align_catalogue()
               if args.save: all_images.write('si-shift')


        #########################################################
        # convolve
        if all_images.suffix_exists('si-conv') and not args.force:
           logging.info('Reuse si-conv images.')
           all_images = lib_fits.AllImages([name.replace('.fits', '-si-conv.fits') for name in all_images.filenames])
        else:
           if args.beam:
              all_images.convolve_to(args.beam, args.circbeam)
           else:
              all_images.convolve_to(circbeam=args.circbeam)
           if args.save: all_images.write('si-conv')
        # regrid
        if all_images.suffix_exists('si-conv-regr') and not args.force:
           logging.info('Reuse si-regr images.')
           all_images = lib_fits.AllImages([name.replace('.fits', '-si-conv-regr.fits') for name in all_images.filenames])
        else:
           all_images.regrid_common(size, args.radec)
           if args.save: all_images.write('si-conv-regr')

        for i, image in enumerate(all_images):
            if args.noise:
               if args.sigma is not None:
                  image.calc_noise(sigma=args.sigma, bg_reg=args.bgreg)  # after mask?/convolution
                  image.blank_noisy(args.sigma)
               else:
                  image.calc_noise() # after mask?/convolution
            if region is not None:
               image.apply_region(region, invert=True) # after convolution to minimise bad pixels
               # all_images.write('blanked')
       
    #########################################################
        # do spdix and write output
        rwcs = pywcs(naxis=2)
        rwcs.wcs.ctype = all_images[0].get_wcs().wcs.ctype
        rwcs.wcs.cdelt = all_images[0].get_wcs().wcs.cdelt
        rwcs.wcs.crval = all_images[0].get_wcs().wcs.crval
        rwcs.wcs.crpix = all_images[0].get_wcs().wcs.crpix
        xsize, ysize = all_images[0].img_data.shape # might be swapped
        regrid_hdr = rwcs.to_header()
        regrid_hdr['NAXIS'] = 2
        regrid_hdr['NAXIS1'] = xsize
        regrid_hdr['NAXIS2'] = ysize
        frequencies = [ image.freq for image in all_images ]
        regrid_hdr['FREQLO'] = np.min(frequencies)
        regrid_hdr['FREQHI'] = np.max(frequencies)
        b = all_images[0].get_beam()
        assert np.all([image.get_beam() == b for image in all_images])
        regrid_hdr['BMAJ'] = b[0]
        regrid_hdr['BMIN'] = b[1]
        regrid_hdr['BPA'] = b[2]
        if args.noise: yerr = np.array([ image.noise for image in all_images ])
        else: yerr = None

        spidx_data = np.nan
        spidx_err_data = np.nan

        flux_output = [image.getflux(args.fluxscaleerr) for image in all_images ]
        flux = [flux_output[i][0] for i in range(len(flux_output)) ]
        flux_err = [flux_output[i][1] for i in range(len(flux_output)) ]
        alpha, error, = twopoint_spidx_bootstrap(frequencies, flux, flux_err, niter=10000)
        alpha = np.log10(flux[0]/flux[1])/np.log10(frequencies[0]/frequencies[1])
        spidx.append(alpha)
        err.append(error)
        print(flux)
        exit()
        
   data = np.column_stack([df.name[:len(spidx)], spidx, err])
   w_file = '/Users/marco/Desktop/The_Master/PhD/GRG Project/File/gmrt-hba_int_spidx.txt'
   spidx_file = open(w_file, "a+")
   np.savetxt(spidx_file, data, fmt='%s')
   spidx_file.close()

 