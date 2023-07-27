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
bootesfile = main + '/BLDF/File/bootes_RGs_catalogue.csv'
elaisfile = main + '/ELDF/File/en1_RGs_catalogue.csv'
lhfile = main + '/LHLDF/File/lh_RGs_catalogue.csv'

bootesdf = pd.read_csv(bootesfile, delimiter=',', usecols=('name', 'z', 'lls'))
elaisdf = pd.read_csv(elaisfile, delimiter=',', usecols=('name', 'z', 'lls'))
lhdf = pd.read_csv(lhfile, delimiter=',', usecols=('name', 'z', 'lls'))

df = pd.concat([bootesdf, elaisdf, lhdf])

df = df[df.lls >= 0.5]
df = df[df.z <= 0.7]
print(df.lls)
exit()

fig = plt.figure(figsize = (6, 6))
fig.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)

x = [0,np.max(df.z)]
y = [0.7, 0.7]

plt.scatter(df.z, df.lls, marker='.', s = 20, edgecolors='black', facecolors='black')
plt.plot(x,y)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('z', fontsize=15)
plt.ylabel('LLS [Mpc]', fontsize=15)
plt.legend(fontsize=12)
plt.yscale('log')
plt.xscale('log')
plt.show()