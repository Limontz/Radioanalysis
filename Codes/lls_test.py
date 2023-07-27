import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
import scipy.stats
import seaborn as sns
from lib_linearfit import linear_fit, linear_fit_bootstrap
from statsmodels.base.model import GenericLikelihoodModel

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project/'
bootes_file = main + 'BLDF/File/bootes_RGs_catalogue.csv'
en1_file = main + 'ELDF/File/en1_RGs_catalogue.csv'
lh_file = main + 'LHLDF/File/lh_RGs_catalogue.csv'
heinz_file = main + 'File/Heinz_RGs_list.csv'
cdf_file = main + 'File/Heinz_cumulative_function.csv'


bootes_df = pd.read_csv(bootes_file, delimiter=',', usecols=['lls'])
en1_df = pd.read_csv(en1_file, delimiter=',', usecols=['lls'])
lh_df = pd.read_csv(lh_file, delimiter=',', usecols=['lls'])

lh_df = pd.read_csv(lh_file, delimiter=',', usecols=['lls'])

df = pd.concat([bootes_df, en1_df])
df = pd.concat([df, lh_df])


df = pd.read_csv(heinz_file, delimiter=',', usecols=['lls'])
cdf = pd.read_csv(cdf_file, delimiter=',')

min = 0.7
max = 1.01
# df = df[df.lls >= min]
# df = df[df.lls < max]

# cdf = cdf[cdf.lls >= min]
# cdf = cdf[cdf.lls < max]

# print(scipy.stats.kstest(cdf.lls, 'uniform', args=(min, max)))
deg_rad_conv = np.pi/180
rad_arcsec_conv = 3600 * 180 / (np.pi)
rad_arcmin_conv = 60 * 180 / (np.pi)

# p = sns.histplot(data=df.lls, bins = 20, zorder = 2, stat='probability')

exp_est = 1
x0_est = 0.5
scale_est = 0.5

def genpareto_fit(x, exp_est, x0_est, scale):
    pl_fit = scipy.stats.genpareto.fit(x, exp_est, loc=x0_est, scale = scale)
    # print(pl_fit)
    x_arr = np.linspace(np.min(x), np.max(x), 1000)
    pdf =  scipy.stats.genpareto.pdf(x_arr, c = (pl_fit[0]), loc = pl_fit[1], scale = pl_fit[2])
    cdf =  scipy.stats.genpareto.cdf(x_arr, c = (pl_fit[0]), loc = pl_fit[1], scale = pl_fit[2])

    test = []
    for i in range(1000):
        r = scipy.stats.genpareto.rvs(size=1000, c = pl_fit[0], loc = pl_fit[1], scale = pl_fit[2] )
        test.append(scipy.stats.kstest(x, r)[1])
    kstest = np.mean(test)
    return(pl_fit,x_arr, pdf, cdf, kstest)

def exp_fit(x, x0_est, scale):
    pl_fit = scipy.stats.expon.fit(x, scale = scale, loc=x0_est)
    # print(pl_fit)
    x_arr = np.linspace(np.min(x), np.max(x), 1000)
    pdf = scipy.stats.expon.pdf(x_arr, loc = pl_fit[0], scale = pl_fit[1])
    cdf = scipy.stats.expon.cdf(x_arr, loc = pl_fit[0], scale = pl_fit[1])

    test = []
    for i in range(1000):
        r = scipy.stats.expon.rvs(size=1000, loc = pl_fit[0], scale = pl_fit[1] )
        test.append(scipy.stats.kstest(x, r)[1])
        # print(test)
    kstest = np.mean(test)
    return(pl_fit,x_arr, pdf, cdf, kstest)

def gumbel_fit(x, x0_est, scale):
    pl_fit = scipy.stats.gumbel_r.fit(x, loc=x0_est, scale = scale)
    # print(pl_fit)
    x_arr = np.linspace(np.min(x), np.max(x), 1000)
    pdf = scipy.stats.gumbel_r.pdf(x_arr, loc = pl_fit[0], scale = pl_fit[1])

    test = []
    for i in range(1000):
        r = scipy.stats.gumbel_r.rvs(size=1000, loc = pl_fit[0], scale = pl_fit[1] )
        test.append(scipy.stats.kstest(x, r))
    kstest = np.mean(test)
    return(pl_fit,x_arr, pdf, kstest)

def number_density(min, fit, coverage, N, zmax):

    volume = (4./3.) * np.pi * cosmo.comoving_distance(zmax)**3
    prob = 1 - scipy.stats.expon.cdf(min, loc = fit[0], scale = fit[1])
    density = N*prob
    # print(density)

    z = np.linspace(0,3, 100)
    number = density*(cosmo.comoving_distance(z)/(0.88*cosmo.comoving_distance(zmax)))**3
    # print(number)
    return(number, z)



pareto_fit, pareto_x_arr, pareto_pdf, pareto_cdf, pareto_kstest = genpareto_fit(df.lls, exp_est, x0_est, scale_est)
exp_fit, exp_x_arr, exp_pdf, exp_cdf, exp_kstest = exp_fit(df.lls, x0_est, scale_est)
# gumbel_fit, gumbel_x_arr, gumbel_pdf, gumbel_kstest = gumbel_fit(df['lls'], x0_est, scale_est)


plt.plot(pareto_x_arr, 1-pareto_cdf, color='red', linestyle='-', zorder=1, linewidth = 3.0, label = 'Pareto')
plt.plot(exp_x_arr, 1-exp_cdf, color='darkorange', linestyle='dotted', zorder=1, linewidth = 3.0, label = 'Exponential')
plt.plot(cdf.lls, cdf.cdf/np.max(cdf.cdf), linestyle="", marker='.', label = 'Heinz CDF')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()

# number, z = number_density(0.7, exp_fit, 0.88, len(df['z2']), np.max(df['z2']))
# plt.plot(z, number)
# plt.show()
#
# exit()

print('PARETO FIT')
print('fit parameters: c=', pareto_fit[0], ' loc=', pareto_fit[1], ' scale=', pareto_fit[2])
print('pareto ks-test=', pareto_kstest)
print('--------------------------------------------------------------------------------------')

print('EXPONENTIAL FIT')
print('fit parameters; loc=', exp_fit[0], ' scale=', exp_fit[1])
print('exp ks-test=', exp_kstest)
print('--------------------------------------------------------------------------------------')


# print('GUMBEL FIT')
# print('fit parameters; loc=', gumbel_fit[0], ' scale=', gumbel_fit[1])
# print('exp ks-test=', gumbel_kstest)
# print('--------------------------------------------------------------------------------------')

# df1 = dfgrg[dfgrg['z'] < 0.5 ]
# df2 = dfgrg.loc[(dfgrg['z'] >= 0.5) & (dfgrg['z'] <= 1)]
# df3 = dfgrg[dfgrg['z'] > 1 ]

# kstest = scipy.stats.kstest(df2['lls'], df3['lls'])
# print(kstest)



########################### 2nd method ###################################
# print('')
# print('--------------------------2nd method-------------------------------')
# print('')

# loc = 0.15 ; scale = 0.1
# exp_params = scipy.stats.expon.fit(df.lls) # params close to but not the same as (shape, loc, scale)
# pareto_params = scipy.stats.genpareto.fit(df.lls)

# class exp(GenericLikelihoodModel):

#     nparams = 2

#     def loglike(self, params):
#         return scipy.stats.expon.logpdf(self.endog, *params).sum()

# class genpareto(GenericLikelihoodModel):

#     nparams = 3

#     def loglike(self, params):
#         return scipy.stats.genpareto.logpdf(self.endog, *params).sum()


# exp_res = exp(df.lls).fit(start_params=exp_params)
# pareto_res = genpareto(df.lls).fit(start_params=pareto_params)
# # res.df_model = len(params)
# print(exp_res.summary())
# print(pareto_res.summary())

# print('')
# print('-------------------------- 3rd method = curve fit -------------------------------')
# print('')
bin = 30
n = plt.hist(df.lls, bins = bin, cumulative = -1,\
            linewidth = 3.0, color = 'navy', density=True)
delta = n[1][1] - n[1][0]
x = [(n[1][i] + n[1][i+1])/2 for i in range(len(n[1])-1)]
# from scipy.optimize import curve_fit

# def exponenial_func(x, mu, b):
#     return np.exp(-(x- mu)/b) / b

# def genpareto_func(x, c, mu, sigma):
#     return (1+c*(x-mu)/sigma)**(-1-1/c) / sigma

# # print('')
# # print('-------------------------- the largest Radio Galaxy -------------------------------')
# # print('')

# # bra = 217.8547029*u.deg
# # bdec = 34.4599830*u.deg

# # c1 = SkyCoord(dfgrg.ra*u.deg, dfgrg.dec*u.deg, frame='icrs')
# # c2 = SkyCoord(bra, bdec, frame='icrs')

# # dfgrg.sep = sep = c1.separation(c2)

# # dfgrg = dfgrg[dfgrg.sep < 2.2*u.deg]


# # zmax = 2
# # Ntot = 73 * 41253
# # print(Ntot)

# # n = 5

# # Ntot = n*(4.*np.pi/3.)*(100)**3
# # print(Ntot)


# # lmax =  exp_fit[0] +  exp_fit[1]*(np.log(1./ exp_fit[1]) - np.log(1./Ntot))






# exp_popt, exp_pcov = curve_fit(exponenial_func, x, n[0], p0=[50, 100])
# pareto_popt, pareto_pcov = curve_fit(genpareto_func, x, n[0], p0=[0.1, 50, 100])


# # pdf = scipy.stats.genpareto.pdf(pareto_x_arr, scale = popt[1], loc=0.00001, c=popt[0])

# print('exponential curve fit')
# print(exp_popt)
# print(exp_pcov)

# print('')

# print('pareto curve fit')
# print(pareto_popt)
# print(pareto_pcov)




plt.clf()
# fig = plt.figure(figsize = (6, 6))
fig1, ax = plt.subplots(figsize = (6, 6))
fig1.subplots_adjust(top=0.990,bottom=0.125,left=0.165,right=0.990,hspace=0.2,wspace=0.2)
p = sns.histplot(data=df.lls, bins = bin, zorder = 2, stat='probability')
plt.plot(pareto_x_arr, delta*pareto_pdf, color='red', linestyle='-', zorder=3, linewidth = 3.0, label = 'Pareto')
plt.plot(exp_x_arr, delta*exp_pdf, color='darkorange', linestyle='dotted', zorder=3, linewidth = 3.0, label = 'Exponential')
# plt.plot(gumbel_x_arr, delta*gumbel_pdf, color='darkviolet', linestyle='--', zorder=3, linewidth = 3.0, label = 'Gumbel')
# plt.grid(color = 'black', zorder=0, linewidth = 1.5)
plt.xlabel(r'$LLS$  [Mpc]',\
           fontsize = 15) #20 with fraction
plt.ylabel(r'Probability', fontsize = 15) #20 with fraction
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
# ax.patch.set_facecolor('#ababab')
# plt.xscale('log')
# ax.patch.set_alpha(0.3)
plt.legend(fontsize=15)
plt.show()

exit()

# pdf2 = (1+pl_fit[0]*((x_arr - pl_fit[1])/pl_fit[2]))**(-1-(1/pl_fit[0]))
# pdf = (1-0.2 6*x_arr)**(-0.26)
# print(pdf)

# Plot CDF of the original
# p = sns.histplot(data=dfgrg['lls'], stat='density')
# plt.plot(x_arr, pdf, color='red', linestyle='--', label='Original')
# plt.plot(x_arr, pdf2, color='red', linestyle='--')
plt.legend()
# plt.show()

n = plt.hist(df['lls'], bins = 20, cumulative = -1,\
            linewidth = 3.0, color = 'navy', density=True)
# plt.clf()
(a, b, sa, sb) = linear_fit(x= n[1][0:-1], y= n[0], yerr=None, tolog=True)
print(a,b)

p = sns.histplot(data=dfgrg['lls'], bins=15, stat='probability')
# x = np.linspace(np.min(dfgrg['lls']), np.max(dfgrg['lls']),100)
# y = 10**(-0.29) * x**(-2.51)
plt.plot(x_arr, pdf, color='violet', linestyle='--', label='Original')
plt.xscale('log')
# plt.plot(x,y, linewidth = 2.0, color = 'red', linestyle="dashed" )
plt.show()
