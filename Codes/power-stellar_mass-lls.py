import numpy as np
import pandas as pd 
from daily_routine import lls
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from astropy import units as u
import sklearn
from sklearn.linear_model import LinearRegression
from scipy.stats import linregress, ks_2samp
from scipy import stats
import statsmodels.api as sm
from scipy.optimize import curve_fit
import pingouin as pg

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
bootesfile = main + '/BLDF/File/bootes_GRGs_catalogue.csv'
bootesbest = main + '/BLDF/File/Crossmatch/full_Bootes_Best_crossmatch.csv'
bootesfluxfile = main + '/BLDF/File/bootes_GRG_flux.csv'

bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
bootesbestdf = pd.read_csv(bootesbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
bootesflux = pd.read_csv(bootesfluxfile, delimiter=',')

bootesdf = pd.merge(bootesdf, bootesbestdf, on=['name'], how='left')
bootesdf = pd.merge(bootesdf, bootesflux, on=['name'], how='left')

elaisfile = main + '/ELDF/File/en1_GRGs_catalogue.csv'
elaisbest = main + '/ELDF/File/Crossmatch/full_en1_Best_crossmatch.csv'
elaisfluxfile = main + '/ELDF/File/en1_grg_flux.csv'

elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
elaisbestdf = pd.read_csv(elaisbest, delimiter=',', usecols=('name', 'AGN_final', 'Mass_cons','SFR_cons'))
elaisflux = pd.read_csv(elaisfluxfile, delimiter=',')

elaisdf = pd.merge(elaisdf, elaisbestdf, on=['name'], how='left')
elaisdf = pd.merge(elaisdf, elaisflux, on=['name'], how='left')

# df = df[df.lls < 0.7]

lhfile = main + '/LHLDF/File/lh_GRGs_catalogue.csv'
lhbest = main + '/LHLDF/File/Crossmatch/full_lh_Best_crossmatch.csv'
lhfluxfile = main + '/LHLDF/File/lh_grg_flux.csv'

lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))
lhbestdf = pd.read_csv(lhbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))
lhfluxdf = pd.read_csv(lhfluxfile, delimiter=',')

lhdf = pd.merge(lhdf, lhbestdf, on=['name'], how='left')
lhdf = pd.merge(lhdf, lhfluxdf, on=['name'], how='left')

df = pd.concat([bootesdf, elaisdf, lhdf])

df = df[df.SFR_cons > -10]

Mpc_to_m = 3.08e22
Jy_to_W = 1.e-26
alpha = 0.7
df['power'] = 4. * np.pi * (1 + df.z)**(alpha - 1) * (df.flux/1e3) *  ((cosmo.luminosity_distance(df.z) * Mpc_to_m / u.Mpc)**2) * Jy_to_W   #units W
df['sSFR'] = 10**(df.SFR_cons)/df.Mass_cons
df['sSFR'] = 10**(df.SFR_cons)/df.Mass_cons


df_herg = df[df['AGN_final']==1]
df_herg = df[df['AGN_final']==1]

df_lerg = df[df['AGN_final']==0]
df_lerg = df[df['AGN_final']==0]




fig = plt.figure(figsize = (6, 6))
fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

plt.scatter(df.lls, df.power/df.Mass_cons, marker='.', s = 100, edgecolors='black', facecolors='blue', label='HERG')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel(r'$Log(SFR [M_{\odot} \ yr^{-1}]$', fontsize=15)
plt.ylabel(r'$ P_{150} [W \ Hz^{-1}] $', fontsize=15)
plt.legend(fontsize=12)
plt.yscale('log')
plt.xscale('log')
# plt.colorbar()

# ax_histx = fig.add_axes(rect_histx, sharex=ax)
# plt.yticks(fontsize = 12)
# plt.minorticks_off()
# # ax_histx.set_yticks([10, 20])
# # ax_histx.hist(z_spec, bins=15, color = 'red')
# plt.ylabel('PDF', fontsize = 12)

# ax_histy = fig.add_axes(rect_histy, sharey=ax)
# plt.xticks(fontsize = 12)
# plt.minorticks_off()
# # ax_histy.set_xticks([10,20])
# plt.xlabel('PDF', fontsize = 12)

# scatter_hist(x = herg_lls, y = herg_SFR, ax = ax, ax_histx = ax_histx, ax_histy = ax_histy)
# ax_histx.hist(lerg_lls, color = 'orange', histtype = 'step', linewidth = 3.0, density='True')
# ax_histy.hist(lerg_SFR, color = 'orange', histtype = 'step', linewidth = 3.0, density='True', orientation='horizontal' )


plt.show()

# from scipy.stats import kstest

# print(kstest(herg_lls, lerg_lls))




