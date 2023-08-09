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
import seaborn as sns 
from matplotlib.ticker import ScalarFormatter


main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'
bootesfile = main + '/BLDF/File/bootes_RGs_catalogue.csv'
elaisfile = main + '/ELDF/File/en1_RGs_catalogue.csv'
lhfile = main + '/LHLDF/File/lh_RGs_catalogue.csv'

bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))

df = pd.concat([bootesdf, elaisdf, lhdf])
df = df[df.lls>0]
print(np.min(df.lls), np.max(df.lls))
exit()

fig, ax = plt.subplots(figsize = (6, 6))
fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

x = [0.7, 0.7]
y = [0, 700]

p = sns.histplot(data=df.lls, bins = 20, zorder=2, color = 'navy')
plt.plot(x,y, marker='', color="red", linewidth=2.5)
plt.xlabel(r'$LLS$  [Mpc]',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'N', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yscale('log')
# plt.xscale('log')
plt.grid(color = 'white', zorder=0, linewidth = 1.5)
ax.patch.set_facecolor('#ababab')
ax.patch.set_alpha(0.3)
ax.yaxis.set_major_formatter(ScalarFormatter())
plt.show()


