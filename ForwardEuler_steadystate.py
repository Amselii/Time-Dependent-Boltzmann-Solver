# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 10:50:24 2025

@author: ameli
"""

import numpy as np
import matplotlib.pyplot as plt

m = 9.11e-28 #g
q = 4.8e-10 #
n_E = 1e8 #1e3-1e10 cm^3
Gamma = 20 #15-30
K_constant = 4*np.pi*q**4/m *n_E*Gamma #7.34e-10 * n_E * Gamma
sigma0 = 6e-7 


E_min = 1e-50 #low
E_max = 4.8e-9 #3keV
nu = 10000
dE = (E_max - E_min)/nu #0.01

E = np.arange(E_min, E_max, dE)
n_i = 1e10
I = 2.17e-11 #13.6 eV
E_avg = 0.6*I #1.6e-11

z_init = np.zeros_like(E)


def S(E): 
    if 4.6e-9< E< 4.8e-9:
        s = 1 / (4.8e-9 - 4.6e-9)#/100
    else:
        s = 0
    
    return s


def f(E, Ej, E_avgi):
    num = 1
    den = (1 + ((E - Ej)/E_avgi)**2)

    return num/den

def g(Ej, Ii, E_avgi):
    num = 1
    den = (1 + ((Ii + Ej)/E_avgi)**2)

    return num/den

def A(E):
    return np.sqrt(2*E/m)




def numerical(z_init, E, n_i):

    z_new = np.copy(z_init)
    z_new[-1] = 0.0
    dE = E[1] - E[0]
    

    for j in range(len(E)-2, 0, -1):

        Ej = E[j]
        

        integral1 = 0.0
        integral2 = 0.0
        for k in range(j+1, len(E)):

            Ep = E[k]


            if Ep > I:
                integral1 += n_i*sigma0*z_new[k]*f(Ep, Ej, E_avg)*dE
                #print("f:", f(Ep, Ej, E_avg))

            if Ep > (I + Ej):

                integral2 += n_i*sigma0*z_new[k]*g(Ej, I, E_avg)*dE
                #print("g", g(Ej, I, E_avg))

        c2 = 0.0
        
        c2 += n_i * sigma0 * E_avg * np.arctan(Ej/E_avg)
        #print(c2)
        
        c3 = ( m / (2 * Ej) ) * K_constant *(1 / dE + 1 / Ej)
        #print(c3)

        LHS =  c2 + c3 
        c4 = ( m / (2 * Ej) ) * K_constant * z_new[j+1] / dE
        #print(c4)
        z_new[j] = (S(Ej) + integral1 + integral2 + c4) / LHS
        #print(integral2)


    return z_new#, c2, c3, c4, integral1



z_new = numerical(z_init, E, n_i)

print(f'z_new: {z_new}')
#print(f'c2: {c2}')
#print(f'c3: {c3}')
#print(f'c4: {c4}')
#print(f'integral1: {integral1}')

#%%
fit = open("./Fortran/data_10000bins_fortran.txt", 'r')
data = np.loadtxt(fit)

#%%
fig, ax = plt.subplots(layout="constrained")

plt.rcParams.update({'font.size': 15})

plt.plot(E*(6.2e8), np.log10(z_new), label = "Python")#np.log10(z_new))#np.log(np.clip(z_new, 0, None)))#
plt.plot(data[:,0]*(6.2e8), np.log10(data[:,1]), label = "Fortran")
plt.xlabel("E [keV]", fontsize = 20)
plt.ylabel("log(z)", fontsize = 20)
plt.xlim(0, 3)
plt.legend()
#plt.figure(constrained_layout=True)
plt.savefig("Fortran_Python_comparison_10000bins.pdf", dpi = 300)

plt.show()

#print(z_new)