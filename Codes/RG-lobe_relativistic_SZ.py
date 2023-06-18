import scipy.special as sc
from scipy.special import gamma, factorial
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math

main = '/Users/marco/Desktop/The_Master/PhD/GRG Project'

#absolute constants
simga_t = 6.65e-24
me = 9.10e-28
c = 3e10
z = 0.00001
kb = 1.4e-16
Tcmb = 2.73
h = 6.626e-27
Io = 270.33 * (Tcmb / 2.755)**3
p1 = 1 #change it according to the energetic found in X-ray observations

#relative constants
p1_2 = 1
p2 = 10
ne = 2.6 * 1e-3 #change it according to X-ray observations
alpha_r = 0.8


#parameters that depend on relative constants
alpha = 3.0 # 2 * alpha_r + 1
A_p1 = (alpha - 1)/ (p1**(1-alpha) - p2**(1-alpha))
A_p1_2 = (alpha - 1)/ (p1_2**(1-alpha) - p2**(1-alpha))

#parameters of the SZ effect

l = 0.15
tau = 2e-6 * (A_p1/A_p1_2) * (ne/1.e-6) * l #l is in Mpc units

def integrand_beta(t, a, b):

    f = t**(a-1) * (1-t)**(b-1)

    return f

def inc_beta(a, b, x):

    beta = quad(integrand_beta, 0, x, args=(a,b))[0]
    # print(beta)
    # print(quad(integrand_beta, 0, x, args=(a,b))[1])
    return beta

energy = (me * c**2 * (alpha-1) / (6 * (p1**(1-alpha)- p2**(1-alpha)))) *\
         ( (inc_beta((alpha-2)/2, (3-alpha)/2,1/(1+p1**2)))\
          - ( inc_beta((alpha-2)/2, (3-alpha)/2,1/(1+p2**2))) )

print('mean energy:',energy*6.242e8/1.e3, 'Gev')

y_nt = energy * tau / (me*c**2)
print( "non-thermal y:", y_nt)


#photon redistribution function

def momentum_spectrum(p, alpha, p1 , p2):

    A = (alpha - 1)/ (p1**(1-alpha) - p2**(1-alpha))
    f = A * p**(-alpha)

    return f


def integrand_beta(t, a, b):

    f = t**(a-1) * (1-t)**(b-1)

    return f

def inc_beta(a, b, x):

    beta = quad(integrand_beta, 0, x, args=(a,b))[0]

    return beta


def prf_part1(p, alpha):


    g = -inc_beta( (1+alpha)/2, -alpha/2, 1/(1+p**2) )\
        -inc_beta( (3+alpha)/2, -(2+alpha)/2, 1/(1+p**2) )*(7+3*alpha)/(3+alpha)\
        -inc_beta( (5+alpha)/2, -(4+alpha)/2, 1/(1+p**2) )*(12+3*alpha)/(5+alpha)

    return g

def prf_part2(t, p, alpha):

    f = p**(-5-alpha)*( (3/(5+alpha) + 2*p**2/(3+alpha)) * (2*np.arcsinh(p) - abs(np.log(t)))  +\
        abs((1-t)/(1+t)) * (  (1+t**2)/((5+alpha)*2*t) + 5/(5+alpha) + 4*p**2/(3+alpha) + 2*p**4/(1+alpha) ) )


    print(f, (1+t**2)/(2*t + (5+alpha)), p**(-5-alpha), (1+t**2)/((5+alpha)*2*t),  t)
    return f

def photon_red_func(s, alpha, p1, p2):

    A = (alpha - 1)/ (p1**(1-alpha) - p2**(1-alpha))
    t = np.exp(s)

    prf = ( prf_part1(max(p2,np.sqrt(t)/2),alpha) - prf_part1(max(p1,np.sqrt(t)/2),alpha) +\
            prf_part2(t, max(p2,np.sqrt(t)/2),alpha) - prf_part2(t, max(p1,np.sqrt(t)/2),alpha) )

    return prf


# print('********************** TEST  TEST  TEST *********************************')
# a = np.exp(-2*np.arcsinh(10))
# b = np.exp(2*np.arcsinh(10))
# print(quad(prf, np.exp(-2*np.arcsinh(10)), np.exp(2*np.arcsinh(10)), args=(10)))
# print()
# A = -3/(32*1e6) * (-np.log(a/b) + a + b -2)
# B= -3/(32*1e6) * (10+8*10**2+4*10**4) * (1 - a - b +a**2/2 + b**2/2)
# C= -3/(32*1e6) * (1/3 - a**2/2 - b**2/2 + a**3/3 + b**3/3)
# D= (3/(8*10**5)) * ((3+3*10**2+10**4)/np.sqrt(1+10**2)) * (b-a+b**2/2-a**2/2)
# E= -(3/(8*10**5)) * ((3+2*10**2)/(2*10)) * 2 * np.arcsinh(10) * (b-a+b**2/2-a**2/2)
# F= (3/(8*10**5)) * ((3+2*10**2)/(2*10)) * ( 5/2 + a*(np.log(a)-1) + b*(np.log(b)-1) + a**2*(-1/4 + np.log(a)/2) + b**2*(-1/4 + np.log(b)/2) )
#
# print(A, B, C, D, E, F)
# print(A+B+C+D+E+F)
#
# print('*************************************************************************')


#*************************************** START PRF **************************************
s = np.arange(-2*np.arcsinh(10),2*np.arcsinh(10), 0.1)

function = np.zeros(len(s))
function2 = np.zeros(len(s))

# for i in range(len(s)):
#
#     # function[i] = photon_red_func(s[i], 2.6, 1, 1e4)
#     function[i] = photon_red_func(s[i], alpha, 1, 1.e4)
#
# # print(function)
# plt.plot(s, function, label='p1=1')
# plt.show()
# exit()
for i in range(len(s)):

    # function[i] = photon_red_func(s[i], 2.6, 1, 1e4)
    function[i] = photon_red_func(s[i], alpha, 1, 10 )
    # function2[i] = photon_red_func(s[i], alpha, 1, p2)

plt.plot(s, (function), label='prf1')
# # plt.plot(s, (function2), label='prf2')
# # plt.yscale('log')
# plt.ylim(-6,0.5)
plt.legend()
plt.show()


exit()


#*************************************** PRF DONE **************************************


def cmb_equation(s,x):

    Tcmb = 2.73
    h = 6.626e-27
    c = 3e10
    kb = 1.4e-16

    if (x/np.exp(s) > 0.01):
       Io = 2 * (kb*Tcmb)**3 / (h*c)**2 * (x*np.exp(-s))**3/(np.exp(x*np.exp(-s))-1.)
       # Io = np.exp(x*np.exp(-s)/(np.exp(-2*np.arcsinh(1.e5))))
    else:
       Io = 2 * (kb*Tcmb) / (c**2) * x*np.exp(-s)*56.8
    return Io

def integrand(s, x, alpha, p1, p2):

    t = np.exp(s)
    # print('t=',t)
    I = photon_red_func(t, alpha, p1, p2) * cmb_equation(s,x) * t

    return I

freq = np.arange(0, 1000, 10) #GHz
J0 = np.zeros(len(freq))
J1 = np.zeros(len(freq))

for i in range(len(freq)):

    J1[i] = quad(integrand, -2*np.arcsinh(p2), 2*np.arcsinh(p2), args=(freq[i]/56.8, alpha, 1, p2))[0] #* (h*c)**2/ (2*(kb*Tcmb)**3)
    J0[i] = cmb_equation(0,freq[i]/56.8) #* (h*c)**2/ (2*(kb*Tcmb)**3)
    #test[i] = integrand(s[i], 100/56.8, alpha, 1, 1e4)
    # print(test[i])
    # print(test2[i])

# plt.plot(freq/56.8, J1, label='i1')
# plt.plot(freq/56.8, J0, label='i0')
# plt.legend()
# plt.show()
# exit()

g = me * c**2 / energy * (J1-J0)*2 * (h*c)**2/ (2*(kb*Tcmb)**3)
deltaI = np.zeros(len(freq))
plt.plot(freq/56.8, g)
plt.show()
exit()
for i in range(len(freq)):

    deltaI[i] = (2*(kb*Tcmb)**3) / (h*c)**2 * y_nt * g[i] * 1e23 #Jy/sr
    # deltaI[i] = y_nt * g[i]

plt.plot(freq, deltaI*1e3/(60*180/np.pi)**2, label='g')
# plt.plot(freq, abs(deltaI), label='g')
# plt.plot(freq, J0, label='g')
# plt.xlim(0,600)
# plt.yscale('log')
plt.legend()

plt.show()
