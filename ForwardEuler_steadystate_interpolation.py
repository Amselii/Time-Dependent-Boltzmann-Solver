# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 08:49:36 2025

@author: ameli
"""

import numpy as np
import matplotlib.pyplot as plt

m = 9.11e-28 #g
q = 4.8e-10 #
n_E = 1e8 #1e3-1e10 cm^3
Gamma = 20 #15-30
K_constant = 4*np.pi*q**4/m *n_E*Gamma
sigma0 = 6e-7 


E_min = 1e-40 #low
E_max = 4.8e-9 #3keV
nu = 5000
dE = (E_max - E_min)/nu #0.01

E = np.arange(E_min, E_max, dE)
n_i = 1e10
I = 2.17e-11 #13.6 eV
E_avg = 0.6*I #1.6e-11

z_init = np.zeros_like(E)


def S(E): 
    if 4.6e-9< E< 4.8e-9:
        s = 1 / (4.8e-9 - 4.6e-9)
    elif 3.0e-9 < E < 4.6e-9:
        s = np.exp(-(4.6e-9 - E)/0.1e-9) / (4.8e-9 - 4.6e-9) #0
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
        
        c2 = 0.0
        
        c2 += n_i * sigma0 * E_avg * np.arctan(Ej/E_avg)
        #print(c2)
        
        c3 = ( m / (2 * Ej) ) * K_constant *(1 / dE + 1 / Ej)
        #print(c3)

        LHS =  c2 + c3 
        c4 = ( m / (2 * Ej) ) * K_constant * z_new[j+1] / dE
        

        integral1 = 0.0
        integral2 = 0.0
        for k in range(j+1, len(E)):

            Ep = E[k]
            


            #if Ep > I:
                #integral1 += n_i*sigma0*z_new[k]*f(Ep, Ej, E_avg)*dE
                
                


            if Ep > (I + Ej):

                integral2 += n_i*sigma0*z_new[k]*g(Ej, I, E_avg)*dE
                #print("g", g(Ej, I, E_avg))

        #z_new[j] = (S(Ej) + integral1 + integral2 + c4) / LHS  
        
        epsilon_r = E[j+1]-Ej
        epsilon_l = 0
       

        integral_ln = n_i*sigma0*(E_avg**2*np.log(epsilon_r**2+E_avg**2)/2 -  E_avg**2*np.log(epsilon_l**2+E_avg**2)/2)
        
        integral_arctan = n_i*sigma0*(E_avg*np.arctan(epsilon_r/E_avg) - E_avg*np.arctan(epsilon_l/E_avg))
        
        #print(integral_ln)#+ integral_arctan)
        
        A_matrix = np.array([[epsilon_l, 1, -1], 
             [epsilon_r, 1, 0], 
             [-integral_ln, -integral_arctan, LHS]])
        
        
        #x = [[A_const, b_const, z_new[j]]]
        B_matrix = np.array([0, z_new[j+1], S(Ej) + integral2 + c4])
        
        x = np.linalg.solve(A_matrix, B_matrix)

        z_new[j] = x[2]

    return z_new



z_new = numerical(z_init, E, n_i)

# print(f'z_new: {z_new}')
#%%

print(z_new)

fit10000 = open("data_10000bins_fortran.txt", 'r')
data10000 = np.loadtxt(fit10000) 

print(data10000)

#%%
plt.plot(data10000[:,0]*(6.2e8), np.log10(data10000[:,1]), label = "10000 bins")
plt.plot(E*(6.2e8), np.log10(z_new))

plt.xlabel("E [keV]")
plt.ylabel("log(z)")
plt.xlim(0, 3)
# plt.savefig("Distribution_steadystate_5000bins.pdf")
plt.show()

#print(z_new)