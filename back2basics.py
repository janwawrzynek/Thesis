#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 14:12:08 2025

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
m_N = 0.1 # HNL Mass in GeV
m_a = 1e-6 # ALP mass in GeV
mixing_parameter_squared = 1e-10 # active-sterile neutrino sector mixing
f_a = 1e3 # ALP decay constant scenario 1
f_M = 0.1 # in GeV

# Integration limits in z corresponding to t=1e-5 s and t=1 s
#z_init = 0.3675    # from t≈1e-5 s (see text)
#z_final = 116.3    # from t=1 s
#z_final = 500
z_init = 1e-4
z_final = 1e2
#Relativistic Degrees of freedom
g_star = 10.75
mass_temp_ratio=  np.linspace(z_init,z_final, 10000)
s_physical =  4.1*m_N**2 # minimum centre of mass energy value required, mandelstam s

#defining HNL mass to temperature ratio

    

def H(z):
    H = (np.pi/M_pl)*(np.sqrt(g_star/90))*((m_N/z)**2)
    return H


    
#entropy density s
def entropy_density(z):
    entropy_density = ((2*(np.pi**2))/(45))*g_star*((m_N/z)**3)
    return entropy_density

'''
def equilibrium_species_density(g,mass,z):
    species_density = g*(mass**3)*(((1/(2*np.pi*z))**1.5))*np.exp(-z)
    return species_density 
#n_N_eq = equilibrium_species_density(2,0.1,z)    
    
'''
def equilibrium_species_density(g,mass,z):
    species_density = g*(mass**3)*(1/(2*(np.pi**2)*(z)))*kn(2,z)
    return species_density 

# equation 2.23, Approximated Total HNL decay to SM
#f_M typical decay constant of a pseudoscalar or vector meson O(0.1GeV)
def decay_width_N_SM(f_M,m_N,mixing_parameter_squared):
    decay_width_N_SM = (((6*(f_M**2))/np.pi)+((m_N**2)/(20*(np.pi**3))))*(m_N**3)*(mixing_parameter_squared)*(fermi_constant**2)
    return decay_width_N_SM
#print(decay_width_N_SM(f_M,m_N,mixing_parameter), "decay")
#now assuming that the particles are relativistic we can obtain the thermally averaged interaction density, as defined in Deppisch 2024, eqn 3.2
#need n_x_eq, ie the HNL equilibrium density

def decay_width_N_av(f_M,m_N,mixing_parameter_squared, m_a):
    decay_width_N_av = mixing_parameter_squared*(m_N**3)*(1/(4*np.pi*(f_a**2)))*(np.sqrt(1+ (m_a/m_N)**2))*((1-(m_a/m_N)**2)**1.5)
    return decay_width_N_av

# This is the part outside the main brackets
def cross_section_aa_NN(m_N,m_a,mixing_parameter, s, f_a):
    
    # Check kinematic thresholds: s must be > max(4m_a², 4m_N²)
    s_threshold = max(4 * m_a**2, 4 * m_N**2)
    if s <= s_threshold:
        return 0.0

    # Precompute common terms
    s_minus_4ma2 = s - 4 * m_a**2
    s_minus_4mN2 = s - 4 * m_N**2
    sqrt_s_minus_4ma2 = np.sqrt(s_minus_4ma2)
    sqrt_s_minus_4mN2 = np.sqrt(s_minus_4mN2)
    
    # Term B: Fraction term
    numerator_B = m_a**4 - 4 * m_a**2 * m_N**2
    denominator_B = numerator_B + m_N**2 * s
    term_B = -1 - numerator_B / denominator_B

    # Term C: coth term
    numerator_C = s - 2 * m_a**2
    denominator_C = sqrt_s_minus_4ma2 * sqrt_s_minus_4mN2
    acoth_arg = numerator_C / denominator_C
    #
    if abs(acoth_arg) <= 1:
        return 0.0

    acoth_val = 0.5 * np.log((acoth_arg + 1) / (acoth_arg - 1))
    term_C = (2 * numerator_C * acoth_val) / denominator_C

    # Combine terms B and C
    bracket_term = term_B + term_C

    # Prefactor
    prefactor = (4 * m_N**4 * sqrt_s_minus_4ma2 ) / (np.pi * f_a**4 * s**1.5)

    # Full cross-section
    cross_section = prefactor * bracket_term
    return cross_section


from scipy.integrate import quad
from scipy.special import kn, k1

def gamma_aa_NN(z, m_N, m_a, mixing_parameter, f_a):

   
    s_min = max((2*m_a)**2, (2*m_N)**2)
    

    s_max = (10*(T))**2  # Safe upper bound where K₁ decays sufficiently
    

    def lambda_func(a, b, c):
        return (a - b - c)**2 - 4*b*c
    

    def reduced_sigma(s):

        sigma_val = cross_section_aa_NN(m_N, m_a, mixing_parameter, s, f_a)
        

        lam = lambda_func(1, m_a**2/s, m_a**2/s)
  
        return 2 * s * sigma_val * lam
    

    def integrand(s):
        # Handle s near zero to avoid numerical issues
        if s <= 0:
            return 0.0
        return np.sqrt(s) * reduced_sigma(s) * k1(np.sqrt(s)/(m_n/z))
    

    # below format ignores error value as integrand produces result and error
    integral, _ = quad(integrand, s_min, s_max, 
                      epsabs=1e-10, epsrel=1e-6, 
                      limit=1000)
    

    gamma = (m_N/z) / (64 * np.pi**4) * integral
    return gamma




#implementing 3.3 of new paper, where K refer to the Bessel functions.
def thermally_averaged_decay_width_density_N_SM(n_N_eq, decay_width_N_SM,z):
    decay_width_N_SM = decay_width_N_SM(f_M,m_N,mixing_parameter_squared)
    thermally_averaged_interaction_density_N_SM = n_N_eq*(kn(1,z)/kn(2,z))*decay_width_N_SM
    return thermally_averaged_interaction_density_N_SM


def dY_N_dz(z, Y_N):
    #dY_N_dz =  -(1 / (z * H_list * s_list)) * gamma_list * ((Y_N / Y_N_eq_list) - 1)
    
    H_ind =  H(z)
    s_ind = entropy_density(z)
    equilibrium_species_density_ind = equilibrium_species_density(2, 0.1, z)
    Y_N_eq_ind = equilibrium_species_density_ind / s_ind
    gamma_ind = thermally_averaged_decay_width_density_N_SM(equilibrium_species_density_ind, decay_width_N_SM, z)
    
    
    dY_N_dz =  -(1 / (z * H_ind * s_ind)) * gamma_ind * ((Y_N / Y_N_eq_ind) - 1)
    return dY_N_dz



H_list = []
s_list = []
n_N_eq_list = []
Y_N_eq_list = []
gamma_list = []
Y_list = []
gamma_n_N_eq_ratio_list =  []
gamma_N_av_n_N_eq_ratio_list = []
gamma_aa_NN_list = []

for z in mass_temp_ratio:
    H_value = H(z)
    s_value = entropy_density(z)
    n_N_eq_value = equilibrium_species_density(2, 0.1, z)
    Y_N_eq_value = n_N_eq_value / entropy_density(z)
    gamma_value = thermally_averaged_decay_width_density_N_SM(n_N_eq_value, decay_width_N_SM, z)
    #Below I am calculating the ratio of 
    gamma_n_N_eq_ratio_value =  decay_width_N_SM(f_M,m_N,mixing_parameter_squared)*(kn(1,z)/kn(2,z))
    gamma_N_av_n_N_eq_ratio_value =  decay_width_N_av(f_M,m_N,mixing_parameter_squared, m_a)*(kn(1,z)/kn(2,z))

     #adding to the lists of values
    H_list.append(H_value)
    s_list.append(s_value)
    n_N_eq_list.append(n_N_eq_value)
    Y_N_eq_list.append(Y_N_eq_value)
    gamma_list.append(gamma_value)
    gamma_n_N_eq_ratio_list.append(gamma_n_N_eq_ratio_value)
    gamma_N_av_n_N_eq_ratio_list.append(gamma_N_av_n_N_eq_ratio_value)
    
#was here


plt.loglog(mass_temp_ratio, gamma_N_av_n_N_eq_ratio_list)
plt.loglog(mass_temp_ratio, gamma_n_N_eq_ratio_list)
plt.loglog(mass_temp_ratio, H_list)
plt.title("gamma_n_N_eq_ratio")
plt.xlabel("z = m_N /T")
plt.ylabel("gamma_n_N_eq_ratio")
plt.show()

#need to fix the definitions of z

    