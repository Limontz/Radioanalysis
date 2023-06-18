#Plot WISE colour-colour plot

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
import math
import pandas as pd

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
magnitude_file = main + '/BLDF/File/Wise_magnitude.grg'
file = main + '/BLDF/File/grgtable2.csv'
flux_file = main + '/BLDF/File/final_bootes_flux.grg'



dfgrg = pd.read_csv(file, delimiter=None, usecols=('name','z','type'))

with open(magnitude_file, 'r') as data:
        tot = []
        for line in data:
            if not line.startswith("#"):
               #p = line.split()
               tot.append(line)

data = np.loadtxt(flux_file, comments='#', delimiter=None, usecols = (1,2))
name = np.loadtxt(flux_file,  dtype='str', comments='#', delimiter=None, usecols = (0))
# dfflux = pd.DataFrame({'name':name, 'flux':data[:,0]})

# df1 = pd.DataFrame({'name': np.concatenate((name,name2))})
# df1.name = 'J' + df1.name


# Bootes part

dfflux = pd.DataFrame({'name':name, 'flux':(data[:,0]), 'fluxerror':(data[:,1])})
dfflux.name = 'J' + dfflux.name

W1 = []
W1_sigma = []
W1_snr = []

W2 = []
W2_sigma = []
W2_snr = []

W3 = []
W3_sigma = []
W3_snr = []

W4 = []
W4_sigma = []
W4_snr = []


CW1 = []
CW1_sigma = []
CW1_snr = []

CW2 = []
CW2_sigma = []
CW2_snr = []

wname = []

for i in range(1, len(tot), 2):
    p = tot[i].split()
    s = tot[i-1].split()

    wname.append(str(s[0]))

    W1.append(float(s[1]))
    W1_sigma.append(float(s[2]))
    W1_snr.append(float(s[3]))

    W2.append(float(s[4]))
    W2_sigma.append(float(s[5]))
    W2_snr.append(float(s[6]))

    W3.append(float(s[7]))
    W3_sigma.append(float(s[8]))
    W3_snr.append(float(s[9]))

    W4.append(float(s[10]))
    W4_sigma.append(float(s[11]))
    W4_snr.append(float(s[12]))


    CW1.append(float(p[0]))
    CW1_sigma.append(float(p[1]))
    CW1_snr.append(float(p[2]))

    CW2.append(float(p[3]))
    CW2_sigma.append(float(p[4]))
    CW2_snr.append(float(p[5]))

#Luminosity estimation at 22 microm
df1 = pd.DataFrame({'name':wname})
df1.name = 'J' +df1.name
dfwise = pd.DataFrame({'name':df1.name, 'W1':W1, 'W1sigma': W1_sigma, 'W1snr':W1_snr, 'W2':W2, 'W2sigma': W2_sigma, 'W2snr':W2_snr, 'W3':W3, 'W3sigma': W3_sigma, 'W3snr':W3_snr, 'W4':W4, 'W4sigma': W4_sigma, 'W4snr':W4_snr})

dfgrg = pd.merge(dfgrg, dfflux, on=['name'], how='left')
# dfgrg = pd.merge(dfgrg, dfflux2, on=['name'], how='left')
dfgrg = pd.merge(dfgrg, dfwise, on=['name'], how='left')
# print(dfgrg)


cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
Mpc_to_cm = 3.08 * 1e24

c = 3e10
wavelength = 22e-3

coeff = 8.363
IR_flux = coeff * 10**(-dfgrg.W4/2.5)
err_IR_flux = coeff * 10**(-dfgrg.W4/2.5) * np.log(10) * dfgrg.W4sigma / 2.5
L = 4. * np.pi * (1 + dfgrg.z)**(-1) * IR_flux * (cosmo.luminosity_distance(dfgrg.z)\
    * Mpc_to_cm / u.Mpc)**2 * 1e-23 #units erg/s/Hz
vL = L * (c/wavelength)
err_vL = 4. * np.pi * (1 + dfgrg.z)**(-1) * (cosmo.luminosity_distance(dfgrg.z)\
         * Mpc_to_cm / u.Mpc)**2 * 1e-23 * err_IR_flux * (c/wavelength)

radio_power = 4. * np.pi * (1 + dfgrg.z)**(- 1) * (cosmo.luminosity_distance(dfgrg.z) * Mpc_to_cm / u.Mpc)**2\
              * (dfgrg.flux / 1.e3) * 1e-23 #units erg/s/HZ
vP = 1.5e8 * radio_power
err_vP =  4. * np.pi * (1 + dfgrg.z)**(- 1) * (cosmo.luminosity_distance(dfgrg.z) * Mpc_to_cm / u.Mpc)**2\
          *(dfgrg.fluxerror / 1.e3) * 1e-23 * 1.5e8


# end of Bootes

#elais part

elaisfile = main + '/ELDF/File/GRG_elais.csv'
elaisfile_flux = main + '/ELDF/File/elais_grg_flux.txt'
elaisfile_wise = main + '/ELDF/File/WISE_grg_crossmatch.csv'

df_elais = pd.read_csv(elaisfile, delimiter=',', usecols=['NAME', 'myra', 'mydec', 'zmean'])
df_elais_wise = pd.read_csv(elaisfile_wise, delimiter=',', usecols=['myra','w1','w1std','w2','w2std','w3','w3std','w4','w4std'])
data = np.loadtxt(elaisfile_flux, comments='#', delimiter=' ', usecols=(1,2))
name = np.loadtxt(elaisfile_flux, dtype='str', comments='#', delimiter=None, usecols = (0))

df_elais_flux = pd.DataFrame({'NAME':name, 'flux':data[:,0], 'fluxerr':data[:,1]})

df_elais = pd.merge(df_elais,df_elais_wise,on='myra', how='left')
df_elais = pd.merge(df_elais,df_elais_flux[['NAME','flux','fluxerr']], on='NAME', how='left')

df_elais = df_elais.rename(columns={"zmean": "z"})

elais_IR_flux = coeff * 10**(-df_elais.w4/2.5)
elais_err_IR_flux = coeff * 10**(-df_elais.w4/2.5) * np.log(10) * df_elais.w4std / 2.5
elais_L = 4. * np.pi * (1 + df_elais.z)**(-1) * elais_IR_flux * (cosmo.luminosity_distance(df_elais.z)\
          * Mpc_to_cm / u.Mpc)**2 * 1e-23 #units erg/s/Hz
elais_vL = elais_L * (c/wavelength)
elais_err_vL = 4. * np.pi * (1 + df_elais.z)**(-1) * (cosmo.luminosity_distance(df_elais.z)\
         * Mpc_to_cm / u.Mpc)**2 * 1e-23 * elais_err_IR_flux * (c/wavelength)

elais_radio_power = 4. * np.pi * (1 + df_elais.z)**(- 1) * (cosmo.luminosity_distance(df_elais.z) * Mpc_to_cm / u.Mpc)**2\
              * (df_elais.flux / 1.e3) * 1e-23 #units erg/s/HZ
elais_vP = 1.5e8 * elais_radio_power
elais_err_vP =  4. * np.pi * (1 + df_elais.z)**(- 1) * (cosmo.luminosity_distance(df_elais.z) * Mpc_to_cm / u.Mpc)**2\
          *(df_elais.fluxerr / 1.e3) * 1e-23 * 1.5e8


# building up the arrays from the 3 fields

IR_flux = np.concatenate((IR_flux, elais_IR_flux), axis=None)
err_IR_flux = np.concatenate((err_IR_flux, elais_err_IR_flux), axis=None)
L = np.concatenate((L, elais_L), axis=None)
vL = np.concatenate((vL, elais_vL),axis=None)
err_vL = np.concatenate((err_vL, elais_err_vL),axis=None)

radio_power = np.concatenate((radio_power, elais_radio_power),axis=None)
vP = np.concatenate((vP, elais_vP),axis=None)
err_vP =  np.concatenate((err_vP, elais_err_vP),axis=None)

w4 = np.concatenate((dfgrg.W4,df_elais.w4),axis=None)
w4std = np.concatenate((dfgrg.W4sigma,df_elais.w4std),axis=None)

plt.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)


err_vL = vL/5
#**************** P_151MHz - P_IR RELATION ***************************

for i in range(len(w4)):

    if (math.isnan(w4std[i]) == True):

        QSO1 = vP[i]
        QSO2 = vL[i]

        # if (dfgrg.type[i] == 'Q'):
        #
        #    QSO1 = vP[i]
        #    QSO2 = vL[i]
        #
        #    plt.errorbar(vP[i]*1.e-7, vL[i]*1.e-7, xerr = err_vP[i]*1.e-7, yerr = err_vL[i]*1.e-7,\
        #                 marker = '.', markersize = '5.0', linestyle = '',\
        #                 uplims = True, color = 'red')

        # else:

        upper1 = vP[i]
        upper2 = vL[i]

        plt.errorbar(vP[i]*1.e-7, vL[i]*1.e-7, xerr = err_vP[i]*1.e-7, yerr = err_vL[i]*1.e-7,\
                     marker = '.', markersize = '5.0', linestyle = '',\
                     uplims = True, color = 'black')

    else:

        # if (dfgrg.type[i] == 'Q'):
        #
        #    plt.errorbar(vP[i]*1.e-7, vL[i]*1.e-7, xerr = err_vP[i]*1.e-7, yerr = err_vL[i]*1.e-7,\
        #                 marker = '.', markersize = '5.0', linestyle = '',\
        #                 color = 'red')
        # else:

        plt.errorbar(vP[i]*1.e-7, vL[i]*1.e-7, xerr = err_vP[i]*1.e-7, yerr = err_vL[i]*1.e-7,\
                 marker = '.', markersize = '1.0', linestyle = '',\
                 color = 'black')


x = [4.e37*1.e-7, 9.e42*1.e-7]
y = [5.e43*1.e-7, 5.e43*1.e-7]
plt.plot(QSO1, QSO2, color = 'red', label = 'QSO')
plt.plot(upper1, upper2, color = 'blue')#, label = 'upper limit')
plt.plot(x, y, color = 'orange', linewidth = 2.0)
plt.xlabel(r'$P_{150MHz}$ [W]', fontsize = 15)
plt.ylabel(r'$P_{22 \rm \mu m }$ [W]', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xscale('log')
plt.yscale('log')
plt.xlim(8.e31, 9.e35)
plt.ylim(3.e34, 1.e38)
# plt.legend(fontsize = 15)
plt.show()

# **************** TEST WISEA AND CWISE W12 MAGNITUDES ****************
#
# plt.clf()
#
# x = [12.5, 18.5]
# plt.figure(figsize = (6,6))
# plt.plot(x,x, color = 'gold')
# plt.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
# counter1 = 0.
# counter2 = 0.
#
# plt.errorbar(dfgrg.W1, CW1, xerr = W1_sigma, yerr = CW1_sigma,\
#              linestyle = '', marker = '.', markersize = 5.0, color = 'black', label = 'W1')
# plt.errorbar(W2, CW2, xerr = W2_sigma, yerr = CW2_sigma,\
#              linestyle = '', marker = '.', markersize = 5.0, color = 'blue', label = 'W2')

# for i in range(len(W1)):
#
#     All_error = W1_sigma[i] + W2_sigma[i]
#     C_error = CW1_sigma[i] + CW2_sigma[i]
#
#     if (CW1_snr[i] <=  W1_snr[i]):
#
#        counter1 += 1
#        if (counter1 == 1):
#           plt.errorbar(W1[i], CW1[i], xerr = W1_sigma[i], yerr = CW1_sigma[i],\
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'red', label = 'W1(multiple match)')
#           plt.errorbar(W2[i], CW2[i], xerr = W2_sigma[i], yerr = CW2_sigma[i],\
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'green', label = 'W2(multiple match)')
#
#        else:
#           plt.errorbar(W1[i], CW1[i], xerr = W1_sigma[i], yerr = CW1_sigma[i],\
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'red')
#           plt.errorbar(W2[i], CW2[i], xerr = W2_sigma[i], yerr = CW2_sigma[i],\
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'green')
#     else:
#
#        counter2 += 1
#
#        if (counter2 == 1):
#           plt.errorbar(W1[i], CW1[i] , xerr = W1_sigma[i], yerr = CW1_sigma[i],\
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'black', label = 'W1')
#           plt.errorbar(W2[i], CW2[i], xerr = W2_sigma[i], yerr = CW2_sigma[i],
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'blue', label = 'W2')
#
#        else:
#
#           plt.errorbar(W1[i], CW1[i], xerr = W1_sigma[i], yerr = CW1_sigma[i],\
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'black')
#           plt.errorbar(W2[i], CW2[i], xerr = W2_sigma[i], yerr = CW2_sigma[i],\
#                       linestyle = '', marker = '.', markersize = 5.0, color = 'blue')


# plt.xlabel(r'W1,W2(AllWISE)', fontsize = 15)
# plt.ylabel(r'W1,W2(CWISE)', fontsize = 15)
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)
# plt.legend(fontsize = 14)
#
# plt.show()

#*******************************************************************************************************

#******************************** WISE COLOUR-COLOUR PLOT **********************************************

w1 = np.concatenate((dfgrg.W1,df_elais.w1),axis=None)
w1std = np.concatenate((dfgrg.W1sigma,df_elais.w1std),axis=None)
w2 = np.concatenate((dfgrg.W2,df_elais.w2),axis=None)
w2std = np.concatenate((dfgrg.W2sigma,df_elais.w2std),axis=None)
w3 = np.concatenate((dfgrg.W3,df_elais.w3),axis=None)
w3std = np.concatenate((dfgrg.W3sigma,df_elais.w3std),axis=None)

plt.clf()
#
for i in range(len(w1)):

    AllW12_error = w1std[i] + w2std[i]
    AllW23_error = w2std[i] + w3std[i]
    AllW24_error = w2std[i] + w4std[i]

    if (math.isnan(w3std[i]) == True):


       if (w1[i]-w2[i] >= 0.7):

          plt.errorbar(w2[i] - w3[i], w1[i] - w2[i], xerr = AllW12_error, yerr = AllW23_error,\
                       xuplims = True, marker = '.', markersize = 5.0, color = 'red')

       elif (vL[i] < 5.e43):

             plt.errorbar(w2[i] - w3[i], w1[i] - w2[i], xerr = AllW12_error, yerr = AllW23_error,\
                        xuplims = True, marker = '.', markersize = 5.0, color = 'blue')

       else:

             plt.errorbar(w2[i] - w3[i], w1[i] - w2[i], xerr = AllW12_error, yerr = AllW23_error,\
                          xuplims = True, marker = '.', markersize = 5.0, color = 'black')



    else:

       if (w1[i]-w2[i] >= 0.7):
           Q1 = w1[i] - w2[i]
           Q2 = w2[i] - w3[i]
           plt.errorbar(w2[i] - w3[i], w1[i] - w2[i], xerr = AllW12_error, yerr = AllW23_error,\
                        marker = '.', markersize = 5.0, color = 'red')

       elif (vL[i] < 5.e43):

           LERG1 = w1[i] - w2[i]
           LERG2 = w2[i] - w3[i]

           plt.errorbar(w2[i] - w3[i], w1[i] - w2[i], xerr = AllW12_error, yerr = AllW23_error,\
                        marker = '.', markersize = 5.0, color = 'blue')

       else:

             HERG1 = w1[i] - w2[i]
             HERG2 = w2[i] - w3[i]

             plt.errorbar(w2[i] - w3[i], w1[i] - w2[i], xerr = AllW12_error, yerr = AllW23_error,\
                          marker = '.', markersize = 5.0, color = 'black')

x1 = [0.7, 6.2]
x2 = [2.5, 2.5]
y1 = [1., 1.]
y2 = [-0.5, 2.0]
plt.plot(Q2, Q1, color = 'red', label = 'QSO')
plt.plot(LERG2, LERG1, color = 'blue', label = 'LERG')
plt.plot(HERG2, HERG1, color = 'black', label = 'HERG')
plt.plot(x1, y1, marker = '', linestyle = '-', color = 'orange')
plt.plot(x2, y2, marker = '', linestyle = '-', color = 'orange')
plt.ylabel(r'[3.4]-[4.6](Vega)', fontsize = 15)
plt.xlabel(r'[4.6]-[12](Vega)', fontsize = 15)
plt.ylim(-0.5, 2.0)
plt.xlim(0.7, 6.2)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.legend(fontsize = 15)
plt.show()
