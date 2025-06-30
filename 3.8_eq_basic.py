#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 13:16:49 2025

@author: janwawrzynek
"""

# implementing equations 3.8 and 3.9
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt
# Physical constants (natural units: GeV, s^-1)
M_pl = 2.435e18          # Planck mass [GeV]
fermi_constant = 1.166e-5           # Fermi constant [GeV^-2]
m_N = 10e-1 # HNL Mass in GeV
mixing_parameter = 10e-10 # active-sterile neutrino sector mixing
f_a = 10e12 # ALP decay constant scenario 1
f_M = 0.1 # in GeV

# Integration limits in z corresponding to t=1e-5 s and t=1 s
z_init = 0.3675    # from t≈1e-5 s (see text)
z_final = 116.3    # from t=1 s
T_init, T_final = m_N/z_init, m_N/z_final
T = np.linspace(T_final,T_init,200)

#defining HNL mass to temperature ratio
def z(m_N,T):
    z = m_N/T 
    return z
z_vals = z(m_N,T)
def H(z):
    H = (m_N**2)/((z**2)*M_pl)
    return H

#Relativistic Degrees of freedom
def rel_deg_freedom():
    rel_deg_freedom = 10.75
    return rel_deg_freedom
    
#entropy density s
def entropy_density(z):
    entropy_density = ((2*(np.pi**2))/(45))*rel_deg_freedom()*((m_N/z)**3)
    return entropy_density


def equilibrium_species_density(g,mass,z):
    species_density = g*(mass**3)*(((1/(2*np.pi*z))**1.5))*np.exp(-z)
    return species_density 
n_N_eq = equilibrium_species_density(2,0.1,z_init)    
    
#def comoving_yield(species_density,entropy_density):
#    comoving_yield = species_density/entropy_density
#    return comoving_yield
#HNL comoving yield at thermal equilibrium Y_N_eq
#Y_N_eq = comoving_yield(n_N_eq,entropy_density)
# above didn't work
Y_N_eq = n_N_eq / entropy_density(z_init)

    
# equation 2.23, Approximated Total HNL decay to SM
#f_M typical decay constant of a pseudoscalar or vector meson O(0.1GeV)
def decay_width_N_SM(f_M,m_N,mixing_parameter):
    decay_width_N_SM = (((6*(f_M**2))/np.pi)+((m_N**2)/(20*(np.pi**3))))*(m_N**3)*((abs(mixing_parameter))**2)*(fermi_constant**2)
    return decay_width_N_SM
print(decay_width_N_SM(f_M,m_N,mixing_parameter), "decay")
#now assuming that the particles are relativistic we can obtain the thermally averaged interaction density, as defined in Deppisch 2024, eqn 3.2
#need n_x_eq, ie the HNL equilibrium density


#implementing 3.3 of new paper, where K refer to the Bessel functions.
def thermally_averaged_decay_width_density_N_SM(n_N_eq, decay_width_N_SM,z):
    decay_width_N_SM = decay_width_N_SM(f_M,m_N,mixing_parameter)
    thermally_averaged_interaction_density_N_SM = n_N_eq*(kn(1,z_vals)/kn(2,z_vals))*decay_width_N_SM
    return thermally_averaged_interaction_density_N_SM
gamma = thermally_averaged_decay_width_density_N_SM(n_N_eq, decay_width_N_SM,z)


H_vals = H(z_vals)
s_vals = entropy_density(z_vals)
n_N_eq_vals = equilibrium_species_density(2, 0.1, z_vals)
gamma_vals = thermally_averaged_decay_width_density_N_SM(n_N_eq_vals, decay_width_N_SM, z_vals)


# Interpolate those arrays for use in ODE
from scipy.interpolate import interp1d

H_interp = interp1d(z_vals, H_vals, kind='cubic', fill_value="extrapolate")
s_interp = interp1d(z_vals, s_vals, kind='cubic', fill_value="extrapolate")
n_N_eq_interp = interp1d(z_vals, n_N_eq_vals, kind='cubic', fill_value="extrapolate")
gamma_interp = interp1d(z_vals, gamma_vals, kind='cubic', fill_value="extrapolate")


def dY_N_dz(z, Y_N):
    H_z = H_interp(z)
    s_z = s_interp(z)
    n_N_eq_z = n_N_eq_interp(z)
    gamma_z = gamma_interp(z)
    return -(1 / (z * H_z * s_z)) * gamma_z * ((Y_N / n_N_eq_z) - 1)

# now solve boltzmann equation using scipy ivp 

sol = solve_ivp(dY_N_dz, [z_init, z_final], [Y_N_eq], dense_output=True)
z_vals = sol.t
Y_vals = sol.y[0]

plt.plot(z_vals,Y_vals)
plt.xlim(max(z_vals), min(z_vals))
plt.xlabel("z = m_N /T")
plt.show()

#need to fix the definitions of z


    