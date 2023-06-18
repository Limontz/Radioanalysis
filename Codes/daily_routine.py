
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
import numpy as np
from astropy import units as u

def z_mean(z, err=False):

    mean=np.nanmean(np.array(z))
    diff = z-mean
    if (err==False): error = np.sqrt(np.nansum(diff**2))

    else: error = np.sqrt(np.nansum(err**2))

    z = [x for x in redshifts if ~np.isnan(np.array(x))]
    if z:
        delta = np.max(diff)
        if diff < 0.1: quality = 3
        if diff > 0.1 and diff < 0.3: quality = 2
        if diff > 0.3 and diff < 0.5: quality = 1
        if diff > 0.5: quality = 0
    else:
        diff = ""

    return(mean,error,quality)

def amin_rad_converter(las):

    rad = las * np.pi / (60 * 180)
    return rad

def rad_amin_converter(las):

    amin = las * 180 / (np.pi * 60)
    return amin


def lls(las, z):

    rad = amin_rad_converter(las)
    lls = cosmo.angular_diameter_distance(z) * rad

    return lls

def radio_power(flux, z, alpha=0.7):

    Mpc_to_m = 3.08e22
    Jy_to_W = 1.e-26

    power = 4. * np.pi * (1 + z)**(alpha - 1) * (cosmo.luminosity_distance(z) * Mpc_to_m / u.Mpc)**2\
            * (flux / 1.e3) * Jy_to_W #units W

    return power
