#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 19:13:26 2025

@author: janwawrzynek
"""


#only including the N to SM decay of 3.8

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt
# Physical constants (natural units: GeV, s^-1)
M_pl = 2.435e18          # Planck mass [GeV]
fermi_constant = 1.166e-5           # Fermi constant [GeV^-2]
m_N = 1e-1 # HNL Mass in GeV
mixing_parameter = 1e-10 # active-sterile neutrino sector mixing
f_a = 1e12 # ALP decay constant scenario 1
f_M = 0.1 # in GeV

# Integration limits in z corresponding to t=1e-5 s and t=1 s
z_init = 0.3675    # from t≈1e-5 s (see text)
z_final = 116.3    # from t=1 s
#T_init, T_final = m_N/z_init, m_N/z_final
#T = np.linspace(T_final,T_init,200)
#Relativistic Degrees of freedom
g_star = 10.75
mass_temp_ratio=  np.linspace(z_init,z_final, 500)

#defining HNL mass to temperature ratio

    

def H(z):
    H = (m_N**2)/((z**2)*M_pl)
    return H


    
#entropy density s
def entropy_density(z):
    entropy_density = ((2*(np.pi**2))/(45))*g_star*((m_N/z)**3)
    return entropy_density


def equilibrium_species_density(g,mass,z):
    species_density = g*(mass**3)*(((1/(2*np.pi*z))**1.5))*np.exp(-z)
    return species_density 
#n_N_eq = equilibrium_species_density(2,0.1,z)    
    

    
# equation 2.23, Approximated Total HNL decay to SM
#f_M typical decay constant of a pseudoscalar or vector meson O(0.1GeV)
def decay_width_N_SM(f_M,m_N,mixing_parameter):
    decay_width_N_SM = (((6*(f_M**2))/np.pi)+((m_N**2)/(20*(np.pi**3))))*(m_N**3)*((abs(mixing_parameter))**2)*(fermi_constant**2)
    return decay_width_N_SM
#print(decay_width_N_SM(f_M,m_N,mixing_parameter), "decay")
#now assuming that the particles are relativistic we can obtain the thermally averaged interaction density, as defined in Deppisch 2024, eqn 3.2
#need n_x_eq, ie the HNL equilibrium density


#implementing 3.3 of new paper, where K refer to the Bessel functions.
def thermally_averaged_decay_width_density_N_SM(n_N_eq, decay_width_N_SM,z):
    decay_width_N_SM = decay_width_N_SM(f_M,m_N,mixing_parameter)
    thermally_averaged_interaction_density_N_SM = n_N_eq*(kn(1,z)/kn(2,z))*decay_width_N_SM
    return thermally_averaged_interaction_density_N_SM


def dY_N_dz(z, Y_N):
    #dY_N_dz =  -(1 / (z * H_list * s_list)) * gamma_list * ((Y_N / Y_N_eq_list) - 1)
    dY_N_dz =  -(1 / (z * H_value * s_value)) * gamma_value * ((Y_N / Y_N_eq_value) - 1)
    return dY_N_dz







H_list = []
s_list = []
n_N_eq_list = []
Y_N_eq_list = []
gamma_list = []
Y_list = []

for z in mass_temp_ratio:
    H_value = H(z)
    s_value = entropy_density(z)
    n_N_eq_value = equilibrium_species_density(2, 0.1, z)
    Y_N_eq_value = n_N_eq_value / entropy_density(z)
    gamma_value = thermally_averaged_decay_width_density_N_SM(n_N_eq_value, decay_width_N_SM, z)

     #adding to the lists of values
    H_list.append(H_value)
    s_list.append(s_value)
    n_N_eq_list.append(n_N_eq_value)
    Y_N_eq_list.append(Y_N_eq_value)
    gamma_list.append(gamma_value)

    
Y_N = solve_ivp(dY_N_dz, [z_init, z_final], [1e-10], method = 'RK45')
z_values = Y_N.t
Y_vals = Y_N.y[0]

     
# now solve boltzmann equation using scipy ivp 

plt.plot(z_values ,Y_vals)

plt.xlabel("z = m_N /T")
plt.ylabel("Comoving Yield Y")
plt.show()

#need to fix the definitions of z

    