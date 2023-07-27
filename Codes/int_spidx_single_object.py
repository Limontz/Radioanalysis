if __name__ == '__main__':

    import os, sys, argparse, logging
    import numpy as np
    from astropy.io import fits as pyfits
    from astropy.wcs import WCS as pywcs
    from astropy.coordinates import match_coordinates_sky
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import pyregion
    from lib_linearfit import linear_fit, linear_fit_bootstrap
    from lib_fits import AllImages
    # https://github.com/astrofrog/reproject
    from reproject import reproject_interp, reproject_exact
    from function import SampleImage
    import itertools
    import pandas as pd
    from cutout2D import cutout
    import matplotlib.pyplot as plt
    import os.path

    reproj = reproject_exact
    logging.root.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description='Make spectral index maps, e.g. spidxmap.py --region ds9.reg --noise --sigma 5 --save *fits')
    # parser.add_argument('images', nargs='+', help='List of images to use for spidx')
    parser.add_argument('--ncpu', dest='ncpu', default=1, type=int, help='Number of cpus to use (default: 1)')
    parser.add_argument('--beam', dest='beam', nargs=3, type=float, help='3 parameters final beam to convolve all images (BMAJ (arcsec), BMIN (arcsec), BPA (deg))')
    # parser.add_argument('--region', dest='region', type=str, help='Ds9 region to restrict analysis')
    parser.add_argument('--noiseregion', dest='noiseregion', type=str, help='Ds9 region to calculate rms noise (default: do not use)')
    parser.add_argument('--noisesigma', dest='noisesigma', default=5, type=int, help='Sigma used in the calc_noise function when no region is specified (default: 5)')
    parser.add_argument('--size', dest='size', nargs=2, type=float, help='Size (ra and dec) of final image in degree (example: 3.5 4.0)')
    # parser.add_argument('--radec', dest='radec', nargs=2, type=float, help='RA/DEC where to center final image in deg (if not given, center on first image - example: 32.3 30.1)')
    parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before calculating spidx (default: false)')
    parser.add_argument('--noise', dest='noise', action='store_true', help='Calculate noise of each image, necessary for the error map (default: false)')
    parser.add_argument('--fluxerr', dest='fluxerr', type=float, help='Fractional flux density error to be added in quadrature (default: 0 - example: 0.05 for 5 per cent)')
    parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate results (default: false)')
    parser.add_argument('--sigma', dest='sigma', type=float, help='Restrict to pixels above this sigma in all images')
    parser.add_argument('--circbeam', dest='circbeam', action='store_true', help='Force final beam to be circular (default: False, use minimum common beam area)')
    parser.add_argument('--bootstrap', dest='bootstrap', action='store_true', help='Use bootstrap to estimate errors (default: use normal X|Y with errors)')
    parser.add_argument('--output', dest='output', default='spidx.fits', type=str, help='Name of output mosaic (default: spidx.fits)')
    # parser.add_argument('--source', dest='source', type=str, help='Name of the source')
    #
    args = parser.parse_args()

    main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
    grgfile = main + 'BLDF/File/bootes_GRGs_catalogue.csv'
    # file2 = main + 'File/apertif_which_fits.txt'
    # main_file = main + 'File/grgtable.csv'

    dfgrg = pd.read_csv(grgfile, usecols=('name', 'ra', 'dec'))
    i =14



    reg = main + '/GMRT_proposal/DS9regions/J' + dfgrg.name[i] + '_rightbridge.reg'
    f1 = main + '/BLDF/File/1424+3436hba.fits'
    f2 = main + '/BLDF/File/1424+3436apt.fits'
    # f3 = main + '/File/spidx/Bootes_LoLLS-conv.fits'

        

        
    images = np.array([f1, f2])


     # check input
    if len(images) < 2:
       logging.error('Requires at lest 2 images.')
       sys.exit(1)

    ########################################################
    # prepare images and convolve+regrid if needed
    print('test-test-test', args.size, reg, (dfgrg.ra[i], dfgrg.dec[i]))
    if np.all([os.path.exists(name.replace('.fits', '-conv-regrid.fits')) for name in images]):
        logging.info('Found convolved+regridded image... restoring.')
        all_images = AllImages([name.replace('.fits', '-conv-regrid.fits') for name in images])
        #print(all_images.filenames)
        regrid_hdr = all_images.images[0].img_hdr
        cd1 = abs(regrid_hdr['CDELT1'])
        cd2 = abs(regrid_hdr['CDELT2'])
        #print(cd1, cd2)
        bmaj = regrid_hdr['BMAJ']
        bmin = regrid_hdr['BMIN']
        if ((cd1-cd2)/cd1)>1.0001 and ((bmaj-bmin)/bmin)>1.0001:
            raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

        gfactor = 2.0*np.sqrt(2.0*np.log(2.0))
        barea = 2.0*np.pi*(bmaj/cd1*bmin/cd2)/(gfactor*gfactor) # in pixels
        print('Beam area=', barea,'pixels')
        #logging.info('Beam area is', barea,'pixels')
        resolution = bmin/gfactor
    else:
        
        all_images = AllImages(images)
        all_images.convolve_to(beam=args.beam, circbeam=args.circbeam)
        if args.save: all_images.write('conv')
        regrid_hdr = all_images.regrid_common(size=args.size, region=reg, action='regrid_header')
        if args.save: all_images.write('conv-regrid')
        regrid_hdr = all_images.images[0].img_hdr
        cd1 = abs(regrid_hdr['CDELT1'])
        cd2 = abs(regrid_hdr['CDELT2'])
        #print(cd1, cd2)
        bmaj = regrid_hdr['BMAJ']
        bmin = regrid_hdr['BMIN']
        if ((cd1-cd2)/cd1)>1.0001 and ((bmaj-bmin)/bmin)>1.0001:
            raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

        gfactor = 2.0*np.sqrt(2.0*np.log(2.0))
        barea = 2.0*np.pi*(bmaj/cd1*bmin/cd2)/(gfactor*gfactor) # in pixels
        print('Beam area=', barea,'pixels')
        # logging.info('Beam area is', barea,'pixels')
        resolution = bmin/gfactor


         # #####################################################
         #find+apply shift w.r.t. lowest noise image
        if args.shift:
            all_images.align_catalogue()
        print('')
        print('CONVOLUTION AND REGRIDDING DONE')
        print('')
        # exit()
        # #########################################################
        # only after regrid+convolve apply mask and find noise
        for image in all_images:
            if args.noise:
                if args.sigma is not None:
                 # usually the sigma used for the blanking is rather low, better to increase it for the calc_noise
                    image.calc_noise(sigma=args.noisesigma, bg_reg=args.noiseregion, force_recalc=True) # after convolution
                    image.blank_noisy(args.sigma)
                else:
                    image.calc_noise(force_recalc=True) # after mask?/convolution

            if reg is not None:
                  image.apply_region(reg, invert=True) # after convolution to minimise bad pixels



    #########################################################
    # do spdix and write output
        frequencies = [ image.get_freq() for image in all_images ]

        if args.noise:
            rmserr = np.array([ image.noise for image in all_images ])
        else: yerr = None
        # print('images noise = ', rmserr)

        flux = np.empty(len(all_images))
        flux= np.array([ np.nansum( image.img_data ) / barea for image in all_images ])

        # print('')
        # print('data', ( image.img_data[~np.isnan(image.img_data)] ))
        # print('')
        for i in range(len(flux)):
            if flux[i] == 0:
                flux[i] = args.sigma*rmserr[i]
                print('')
                print('WARNING: Flux = 0 at freq = ', frequencies[i])
                print('')

        apertif_limit = '-'
        lba_limit = '-'

        if (flux[1] == 0 or flux[1] < 4.e-4):
            apertif_limit = 'U'
            flux[1] = args.sigma*4.e-4

        # if (flux[2] < flux[1] or flux[2] == 0):
        #    lba_limit = 'NO'
        #    print('')
        #    print('WARNING: Flux = 0 at freq = ', frequencies[2])
        #    print('')


        #if np.isnan(flux).any() or (np.array(flux) < 0).any(): continue
        # add flux error
        if args.fluxerr:
            yerr = np.sqrt((args.fluxerr*flux)**2+rmserr**2)
        else:
            yerr = rmserr
        #if (np.array(flux) <= 0).any(): continue
        if args.bootstrap:
            (a, b, sa, sb) = linear_fit_bootstrap(x=frequencies, y=flux, yerr=yerr, tolog=True)
        else:
            (a, b, sa, sb) = linear_fit(x=frequencies, y=flux, yerr=yerr, tolog=True)


        flux_perm = list(itertools.combinations(flux, 2))
        freq_perm = list(itertools.combinations(frequencies, 2))
        rmserr_perm = list(itertools.combinations(rmserr, 2))
        # print(flux_perm, rmserr_perm)
        #
        int_a = np.empty(len(flux_perm))
        int_sa = np.empty(len(flux_perm))
        # print(freq_perm)
        for i in range(len(flux_perm)):
            yerr = []
        if args.fluxerr:
            # print(i, ' ',flux_perm[i][0], ' ', rmserr_perm[i][0])
            # print(i, ' ',flux_perm[i][1], ' ', rmserr_perm[i][1])
            yerr_1 = np.sqrt((args.fluxerr*flux_perm[i][0])**2+rmserr_perm[i][0]**2)
            yerr.append(yerr_1)
            yerr_2 = np.sqrt((args.fluxerr*flux_perm[i][1])**2+rmserr_perm[i][1]**2)
            yerr.append(yerr_2)
        else:
            yerr = rmserr[i]
        #if (np.array(flux) <= 0).any(): continue
        if args.bootstrap:
            (int_a[i], b, int_sa[i], sb) = linear_fit_bootstrap(x=freq_perm[i], y=flux_perm[i], yerr=yerr, tolog=True)
        else:
            (int_a[i], b, int_sa[i], sb) = linear_fit(x=freq_perm[i], y=flux_perm[i], yerr=yerr, tolog=True)
        print('')
        print('Integrated spectral index over all images = ', '%.3f'%(a), ' +- ', '%.3f'%(sa))
        print('')
        print('frequency range = ', '%.3f'%(freq_perm[0][0]/1.e6), '-', '%.3f'%(freq_perm[0][1]/1.e6), 'MHz ',\
        '%.3f'%(freq_perm[1][0]/1.e6), '-', '%.3f'%(freq_perm[1][1]/1.e6), 'MHz ',\
        '%.3f'%(freq_perm[2][0]/1.e6), '-', '%.3f'%(freq_perm[2][1]/1.e6), 'MHz ')
        print('spectral index =  ', '%.3f'%(int_a[0]),' +- ', '%.3f'%(int_sa[0]),\
        '     ', '%.3f'%(int_a[1]), ' +- ', '%.3f'%(int_sa[1]), '    ',\
        '%.3f'%(int_a[2]), ' +- ', '%.3f'%(int_sa[2]))


        # os.system('rm *-conv-regrid.fits')
        # os.system('rm lba.fits')
        # os.system('rm apertif.fits')
