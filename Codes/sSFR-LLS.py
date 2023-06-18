import numpy as np
import pandas as pd 
from daily_routine import lls
import matplotlib.pyplot as plt

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
elais_file = main + 'ELDF/File/en1_full_list_RGs.csv'
bootes_file = main + 'BLDF/File/bootes_full_list_RGs.las'
Best_elais_file = main + 'ELDF/File/crossmatch/Best_total_crossmatch.csv'
Best_bootes_file = main + 'BLDF/File/crossmatch/Best_total_crosmmatch.csv'

df_elais = pd.read_csv(elais_file, delimiter=',')

data = np.loadtxt(bootes_file, comments='#', delimiter=None, usecols=(4))
name = np.loadtxt(bootes_file, dtype='str', comments='#', delimiter=None, usecols = (0))

df_bootes = pd.DataFrame({'NAME':name, 'las':data})

df_Best_elais = pd.read_csv(Best_elais_file, delimiter=',')
df_Best_bootes = pd.read_csv(Best_bootes_file, delimiter=',')

df_elais = pd.merge(df_Best_elais,df_elais[['NAME','lls']],on='NAME', how='left')
df_bootes = pd.merge(df_Best_bootes,df_bootes[['NAME','las']],on='NAME', how='left')

df_elais['sSFR'] = 10**(df_elais.SFR_cons)/10**(df_elais.Mass_cons)
df_bootes['sSFR'] = 10**(df_bootes.SFR_cons)/10**(df_bootes.Mass_cons)


df_elais = df_elais[df_elais['SFR_cons'] > -10]
df_bootes = df_bootes[df_bootes['SFR_cons'] > -10]

df_bootes['lls'] = lls(df_bootes['las'], df_bootes['z_best'])

df_herg_elais = df_elais[df_elais['AGN_final']==1]
df_herg_bootes = df_bootes[df_bootes['AGN_final']==1]

df_lerg_elais = df_elais[df_elais['AGN_final']==0]
df_lerg_bootes = df_bootes[df_bootes['AGN_final']==0]

herg_lls = np.concatenate((df_herg_elais.lls, df_herg_bootes.lls),axis=None)
lerg_lls = np.concatenate((df_lerg_elais.lls, df_lerg_bootes.lls),axis=None)

herg_SFR = np.concatenate((df_herg_elais.SFR_cons, df_herg_bootes.SFR_cons),axis=None)
lerg_SFR = np.concatenate((df_lerg_elais.SFR_cons, df_lerg_bootes.SFR_cons),axis=None)

herg_sSFR = np.concatenate((df_herg_elais.sSFR, df_herg_bootes.sSFR),axis=None)
lerg_sSFR = np.concatenate((df_lerg_elais.sSFR, df_lerg_bootes.sSFR),axis=None)



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

plt.scatter(herg_lls, np.log10(herg_sSFR), marker='.', s = 100, edgecolors='black', facecolors='blue', label='HERG')
plt.scatter(lerg_lls, np.log10(lerg_sSFR), marker='.', s = 100, edgecolors='black', facecolors='orange', label='LERG')

plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('LLS [Mpc]', fontsize=15)
plt.ylabel(r'$Log(sSFR [yr^{-1}])$)', fontsize=15)
plt.legend(fontsize=12)
# plt.xscale('log')
# plt.yscale('log')

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

scatter_hist(x = herg_lls, y = np.log10(herg_sSFR), ax = ax, ax_histx = ax_histx, ax_histy = ax_histy)
ax_histx.hist(lerg_lls, color = 'orange', histtype = 'step', linewidth = 3.0, density='True')
ax_histy.hist(np.log10(lerg_sSFR), color = 'orange', histtype = 'step', linewidth = 3.0, density='True', orientation='horizontal' )


plt.show()

from scipy.stats import kstest

print(kstest(herg_sSFR, lerg_sSFR))








