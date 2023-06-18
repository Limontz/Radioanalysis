import numpy as np
import pandas as pd 
from daily_routine import lls
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd 
from daily_routine import lls
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from astropy import units as u
from scipy.stats import linregress, ks_2samp


main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
bootesfile = main + '/BLDF/File/bootes_RGs_catalogue.csv'
bootesbest = main + '/BLDF/File/Crossmatch/full_Bootes_Best_crossmatch.csv'

bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
bootesbestdf = pd.read_csv(bootesbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))

bootesdf = pd.merge(bootesdf, bootesbestdf, on=['name'], how='left')

elaisfile = main + '/ELDF/File/en1_RGs_catalogue.csv'
elaisbest = main + '/ELDF/File/Crossmatch/full_en1_Best_crossmatch.csv'

elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
elaisbestdf = pd.read_csv(elaisbest, delimiter=',', usecols=('name', 'AGN_final', 'Mass_cons','SFR_cons'))

elaisdf = pd.merge(elaisdf, elaisbestdf, on=['name'], how='left')

# df = df[df.lls < 0.7]

lhfile = main + '/LHLDF/File/lh_RGs_catalogue.csv'
lhbest = main + '/LHLDF/File/Crossmatch/full_lh_Best_crossmatch.csv'

lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))
lhbestdf = pd.read_csv(lhbest, delimiter=',', usecols=('name','AGN_final', 'Mass_cons','SFR_cons'))

lhdf = pd.merge(lhdf, lhbestdf, on=['name'], how='left')


df = pd.concat([bootesdf, elaisdf, lhdf])

df = df[df.SFR_cons > -10]

df_herg = df[df['AGN_final']==1]
df_herg = df[df['AGN_final']==1]

df_lerg = df[df['AGN_final']==0]
df_lerg = df[df['AGN_final']==0]


def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    # ax.scatter(x, y, color = 'black', s = 10.0)

    # now determine nice limits by hand:
    # binwidth = 0.25
    # xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    # lim = (int(xymax/binwidth) + 1) * binwidth
    #
    # bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, color = 'navy', density='True', stacked='True')
    ax_histy.hist(y, orientation='horizontal', color = 'navy', density='True')

left, width = 0.165, 0.62
bottom, height = 0.1, 0.68
spacing = 0.005


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes(rect_scatter)

plt.scatter(df_herg.lls, df_herg.SFR_cons, marker='.', s = 100, edgecolors='black', facecolors='blue', label='HERG')
plt.scatter(df_lerg.lls, df_lerg.SFR_cons, marker='.', s = 100, edgecolors='black', facecolors='orange', label='LERG')

plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('LLS [Mpc]', fontsize=15)
plt.ylabel(r'$Log(SFR [M_{\odot} \ yr^{-1}]$)', fontsize=15)
plt.legend(fontsize=12)
plt.xscale('log')

ax_histx = fig.add_axes(rect_histx, sharex=ax)
plt.yticks(fontsize = 12)
plt.minorticks_off()
# ax_histx.set_yticks([10, 20])
# ax_histx.hist(z_spec, bins=15, color = 'red')
plt.ylabel('PDF', fontsize = 12)

ax_histy = fig.add_axes(rect_histy, sharey=ax)
plt.xticks(fontsize = 12)
plt.minorticks_off()
# ax_histy.set_xticks([10,20])
plt.xlabel('PDF', fontsize = 12)

scatter_hist(x = df_herg.lls, y = df_herg.SFR_cons, ax = ax, ax_histx = ax_histx, ax_histy = ax_histy)
ax_histx.hist(df_lerg.lls, color = 'orange', histtype = 'step', linewidth = 3.0, density='True')
ax_histy.hist(df_lerg.SFR_cons, color = 'orange', histtype = 'step', linewidth = 3.0, density='True', orientation='horizontal' )


plt.show()

from scipy.stats import kstest

print(kstest(df_herg.lls, df_lerg.lls))




