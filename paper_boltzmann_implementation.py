#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: janwawrzynek
"""

#Attempt to recreate the papers boltzmann equation solutions in Python

import numpy as np
import matplotlib as plt
import scipy


def s(T):
    s = (2/45)*(np.pi**2)*(g)*(T**3)
    return s 

#gstar is the number of relativistic degrees of freedom. 
#gx is the internal degrees of freedom
# T is the temperature of the thermal bath of SM particles

#comoving yield for any species x

Y = nx/s

#collision operator which takes account of all energy and number changing processes
#def collision_operator(f):
#    C = 
#    return C
    

#Boltzmann equation, Energy density of various particle species due to the expansion of the universe
#and the HNL and axion interaction processes

def energy_density_boltzmann(d_rho_x_dt,H, rho_x,px):
    energy_density_transfer_rate = d_rho_x_dt+3*H(rho_x+px)
    return energy_density_transfer_rate


# boltzmann equation defining number density of various particle species due to expansion of the universe
# and the HNL and axion interaction processes

def number_density_boltzmann(d_n_x_dt,H, n_x):
    number_density_transfer_rate = d_n_x_dt+ 3*H*n_x
    return number_density_transfer_rate


#Transfer rates are obtained from calculating thermally averaged cross-sections and decay rates for each process that species x is involved

def thermally_averaged_interaction_density(cross_section, n_x_eq, decay_rate):
    thermally_averaged_interaction_density = (cross_section)*(n_x_eq**2)
    thermally_averaged_interaction_density = (decay_rate)*(n_x_eq)
    

    
def z(m_x,T):
    z = m_x/T
    return z

def bessel_1(z):
    return K_1

def bessel_2(z):
    return K_2


def thermally_averaged_decay_width_density(n_x_eq,K_1,K_2, zero_temperature_decay_width_x)
    thermally_averaged_decay_width_density = n_x_eq*(K_1/K_2)*zero_temperature_decay_width_x
    return thermally_averaged_decay_width_density

a = 1  

def b(m_x1,s):
    b = ((m_x1)**2)/s
    return b

def c(m_x2,s):
    c = (m_x2**2)/s
    return c



def lamda(a,b,c):
    lamda = (a-b-c)**2  -4*b*c
    return lamda

def reduced_cross_section(s,sigma_s,lamda ):
    reduced_cross_section = 2*s*sigma_s*lamda
    return reduced_cross_section

def s_min(m_x1,m_x2,m_y1,m_y2):
    max([(m_x1+m_x2)**2, (m_y1 + m_y2)**2])
    return s_min

def thermally_averaged_cross_section_density(T,s, reduced_cross_section,K_1,ds):
    #fix integration boundaries 
    thermally_averaged_cross_section_density = (T/(64*(np.pi**2))*scipy.int((s**0.5)*reduced_cross_section*K_1*((s**0.5)/T)*ds)
        
                                                

# Abundances before BBN
# boltzmann equation for HNL number density  before neutrino decoupling

def before_BBN_number_N_boltzmann(decay_rate_N_SM, n_N, n_N_eq, decay_rate_N_av, n_a, n_a_eq, H):
    d_n_N_dt = -decay_rate_N_SM*(n_N-n_N_eq)-decay_rate_N_av*(n_N-n_N_eq*(n_a/n_a_eq))-3*H_n_N
    return d_n_N_dt


# boltzmann equation for Alp number density before neutrino decoupling.

def before_BBN_number_a_boltzmann(decay_rate_a_vv, n_a, n_a_eq, decay_rate_N_av, n_N, n_N_eq, H):
    d_n_a_dt = -decay_rate_a_vv*(n_a-n_a_eq) + decay_rate_N_av*(n_N-n_N_eq*(n_a/n_a_eq)) - 3*H*n_a



    