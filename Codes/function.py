#Function for the grg Project

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
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

from astropy.wcs import WCS as pywcs
from astropy.io import fits
import numpy as np
import os, sys, logging, re
import pandas as pd
import pyregion
#from lib_fits import Image


def flatten(filename, channel=0, freqaxis=0):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    f = fits.open(filename)#, ignore_missing_end=True)

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        pass
        #return f[0].header,f[0].data

    w = pywcs(f[0].header)
    wn = pywcs(naxis=2)

    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]

    header = wn.to_header()
    header["NAXIS"]=2
    header["NAXIS1"]=f[0].header['NAXIS1']
    header["NAXIS2"]=f[0].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    dataslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dataslice.append(np.s_[:],)
        elif i==freqaxis:
            dataslice.append(channel)
        else:
            dataslice.append(0)

    # add freq
    header["FREQ"] = find_freq(f[0].header)

    # add beam if present
    try:
        header["BMAJ"]=f[0].header['BMAJ']
        header["BMIN"]=f[0].header['BMIN']
        header["BPA"]=f[0].header['BPA']
    except:
        pass

    # slice=(0,)*(naxis-2)+(np.s_[:],)*2
    # print(header)
    return header, f[0].data[tuple(dataslice)]

def find_freq(header):
    """
    Find frequency value in most common places of a fits header
    """
    if not header.get('RESTFRQ') is None and not header.get('RESTFRQ') == 0:
        return header.get('RESTFRQ')
    elif not header.get('FREQ') is None and not header.get('FREQ') == 0:
        return header.get('FREQ')
    else:
        for i in range(5):
            type_s = header.get('CTYPE%i' % i)
            if type_s is not None and type_s[0:4] == 'FREQ':
                return header.get('CRVAL%i' % i)

    return None # no freq information found

def correct_beam_header(header):
     """
     Find the primary beam headers following AIPS convenction
     """
     if ('BMAJ' in header) and ('BMIN' in header) and ('PA' in header): return header
     elif 'HISTORY' in header:
         for hist in header['HISTORY']:
             if 'AIPS   CLEAN BMAJ' in hist:
                 # remove every letter from the string
                 bmaj, bmin, pa = re.sub(' +', ' ', re.sub('[A-Z ]*=','',hist)).strip().split(' ')
                 header['BMAJ'] = float(bmaj)
                 header['BMIN'] = float(bmin)
                 header['BPA'] = float(pa)
     return header


def flatten2(filename, channel=0, freqaxis=0):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    f = fits.open(filename)#, ignore_missing_end=True)

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        pass
        #return f[0].header,f[0].data

    w = pywcs(f[0].header)
    wn = pywcs(naxis=2)

    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]

    header = wn.to_header()
    header["NAXIS"]=2
    header["NAXIS1"]=f[0].header['NAXIS1']
    header["NAXIS2"]=f[0].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    dataslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dataslice.append(np.s_[:],)
        elif i==freqaxis:
            dataslice.append(channel)
        else:
            dataslice.append(0)

    # add freq
    header["FREQ"] = find_freq(f[0].header)

    # add beam if present
    try:
        header["BMAJ"]=f[0].header['BMAJ']
        header["BMIN"]=f[0].header['BMIN']
        header["BPA"]=f[0].header['BPA']
    except:
        pass

    # slice=(0,)*(naxis-2)+(np.s_[:],)*2
    # print(header)
    return header, f[0].data[tuple(dataslice)], w


def calculate_flux(imagename, regionfile, nsigma=3.):
    im = SampleImage(imagename)
    print('Calculating flux of ', str(regionfile) )
    im.calc_noise()
    im.set_region(regionfile, individual=True)
    fluxes, errors, upper = im.get_flux(nsigma) #Calculate flux in >=3 sigma pixels only
    return fluxes, errors, upper


class SampleImage(object):

    def __init__(self, imagefile):
        """
        imagefile: name of the fits file
        """
        self.imagefile = imagefile
        self.img_hdr, self.img_data = flatten(self.imagefile)
        self.hdu = fits.PrimaryHDU(self.img_data, self.img_hdr)

        # calculate beam area in pixels
        cd1 = abs(self.img_hdr['CDELT1'])
        cd2 = abs(self.img_hdr['CDELT2'])
        bmaj = self.img_hdr['BMAJ']
        bmin = self.img_hdr['BMIN']

        if ((cd1-cd2)/cd1)>1.0001 and ((bmaj-bmin)/bmin)>1.0001:
                raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

        gfactor = 2.0*np.sqrt(2.0*np.log(2.0))
        self.barea = 2.0*np.pi*(bmaj/cd1*bmin/cd2)/(gfactor*gfactor) # in pixels
        #print('Beam area=',self.barea,'pixels')
        logging.info('Beam area is',self.barea,'pixels')
        self.resolution = bmin/gfactor
        self.mask_noise = None

    def calc_noise(self, niter=1000, eps=None, sigma=5):
        """
        Return the rms of all the pixels in an image
        niter : robust rms estimation
        eps : convergency criterion, if None is 1% of initial rms
        """
        if eps == None: eps = np.nanstd(self.img_data)*1e-3
        data = self.img_data[ ~np.isnan(self.img_data) ] # remove nans
        if len(data) == 0: return 0
        oldrms = 1.
        for i in range(niter):
            rms = np.nanstd(data)
            if np.abs(oldrms-rms)/rms < eps:
                self.noise = rms
                print('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))
                logging.debug('%s: Noise: %.3f mJy/b' % (self.imagefile, self.noise*1e3))
                return rms

            data = data[np.abs(data)<sigma*rms]
            oldrms = rms
        print('Noise=',self.noise)
        exit()
        raise Exception('Noise estimation failed to converge.')



    def set_region(self, regionfile, individual=False):
        self.masks = []
        region = pyregion.open(regionfile).as_imagecoord(self.img_hdr)
        if individual:
            for region_split in region:
                self.masks.append( pyregion.ShapeList([region_split]).\
                get_mask(hdu=self.hdu,shape=np.shape(self.img_data)) )
        else:
            self.masks.append( region.get_mask(hdu=self.hdu,shape=np.shape(self.img_data)) )
        
        return self.masks


    def get_flux(self, nsigma = 3, with_upper_limits =False, upper_limit_sigma = 3):
        """
        nsigma: use only pixels above this sigma
        with_upper_limits: if no detection, set the value at upper_limit_sigma sigma. It also returns a bool array with True for limits
        upper_limit_sigma: numer of sigmas to consider a flux a limit (default: 3)
        """

        # set self.noise
        if self.mask_noise is None:
            noise = self.calc_noise()
        else:
            self.noise = np.nanstd( self.img_data[self.mask_noise] )


        fluxes = []
        errors = []
        for mask in self.masks:
            print('beam area ',self.barea)
            mask = np.logical_and(mask, ~np.isnan(self.img_data))
            #print('non zero',np.count_nonzero(mask))
            p = self.img_data[mask]
            area = len(p)
            print('region area ',area)
            p = p[p > nsigma*self.noise]
            flux = np.nansum(p) / self.barea
            print(flux, (area/self.barea))
            flux = flux * (area/self.barea)
            print(flux)
            # flux = flux * ((self.barea)/area)
            #print(np.shape(self.img_data),np.shape(self.img_data[mask]))
            noise_error = self.noise * np.sqrt( np.count_nonzero(mask) / self.barea )
            print('Assuming 10% calibration error.')
            error = np.sqrt(noise_error**2 + (flux*0.1)**2)
            fluxes.append(flux)
            errors.append(error)

        fluxes = np.array(fluxes)
        errors = np.array(errors)

        if with_upper_limits:
            upper_limits = np.zeros_like(fluxes, dtype=bool)
            is_limit = np.where(fluxes < (upper_limit_sigma * errors))
            upper_limits[ is_limit ] = True
            fluxes[ is_limit ] = upper_limit_sigma*errors[ is_limit ]
            return fluxes, errors, upper_limits
        else:
            return fluxes, errors, False

    def signal_to_noise(self, regionfile, sigma):
        data= self.img_data[ ~np.isnan(self.img_data) ]

        region = self.set_region(regionfile)
        mask = np.logical_and(region, region)
        print(mask[0,:,:].shape, mask[mask==True])
        data = self.img_data[mask[0,:,:]]

        S_N = np.mean(data) / ( (bmaj*bmin) * self.noise * np.sqrt( np.count_nonzero(mask) / self.barea ) )

        print(S_N)
