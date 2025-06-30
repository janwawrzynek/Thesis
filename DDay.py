

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 21:23:51 2025

@author: janwawrzynek
"""

#seed vault changing the relativistic degrees of freedom to be updated inside the z values loop


#only including the N to SM decay of 3.8
#early universe the g_star ratio is around 100, this could explain the issues 

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import jv, yv
import matplotlib.pyplot as plt
import mpmath
# Physical constants (natural units: GeV, s^-1)
M_pl = 2.435e18          # Planck mass [GeV]
fermi_constant = 1.166e-5           # Fermi constant [GeV^-2]
m_N = 0.1 # HNL Mass in GeV
m_a = 1e-6 # ALP mass in GeV
mixing_parameter_squared = 1e-10 # active-sterile neutrino sector mixing
m_v =  1e-9
f_a = 1e3 # ALP decay constant scenario 1
#f_M = 0.1 # in GeV

#meson vector and pseudoscalar masses
#mdxon values for HNL Decays all in GeV
mass_pion_charged = 0.13957 #GeV
mass_pion_neutral = 0.1349766 #GeV
mass_kaon_charged = 0.493677
mass_kaon_neutral = 0.497614 
mass_eta = 0.547862
mass_eta_prime = 0.95778
mass_D_charged = 1.86961
mass_D_strange_charged = 1.96830
mass_omega = 0.78265
mass_phi = 1.019461
mass_rho_charged = 0.77511
mass_rho_neutral = 0.77526
mass_kaon_star_neutral = 0.89581
mass_kaon_star_charged = 0.89166
#meson decay constants
#pseudoscalar mesons for now assume that the neutral and charged mesons have the same decay constants
#all in GeV
decay_constant_pion = 0.130 #GeV
decay_constant_kaon = 0.156 #GeV
decay_constant_D = 0.212 #GeV
decay_constant_D_strange = 0.249 #GeV
#Vector Mesons in GeV^2
decay_constant_rho = 0.171 # GeV^2
decay_constant_omega = 0.155# GeV^2
decay_constant_phi = 0.232# GeV^2
decay_constant_kaon_star = 0.178# GeV^2


#Grouping the pseudoscalars into lists
pseudo_neutral_masses = [mass_pion_neutral, mass_kaon_neutral, mass_eta, mass_eta_prime]


    



# Integration limits in z corresponding to t=1e-5 s and t=1 s
#z_init = 0.3675    # from t≈1e-5 s (see text)
#z_final = 116.3    # from t=1 s
#z_final = 500
z_init = 1e-4
z_final = 1e2
#Relativistic Degrees of freedom
#g_star = 10.75
#g_star = 111

mass_temp_ratio=  np.linspace(z_init,z_final, 100)
s_physical =  4.1*m_N**2 # minimum centre of mass energy value required, mandelstam s

#defining HNL mass to temperature ratio
s_min = max((2*m_a)**2, (2*m_N)**2)
s_max = (10*(m_N/z_final))**2  # Safe upper bound where K₁ decays sufficiently
s_range = np.linspace(s_min, s_max, 50)
m_l = [0.000511, 0.105658, 1.777] #Electron, muon, tau masses in GeV
    

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

#Full Version of N-> SM
sin2w = 0.23126
sinw2 = sin2w
cosw2 = 1 - sin2w

def L_func(x):
    """Eq. (2.13)"""
    return np.log((1 - 3*x**2 - (1-x**2) * (np.sqrt(1 - 4*x**2))) / (x**2 * (1 + np.sqrt(1 - 4*x**2))))

def f1(x):
    """Eq. (2.12), first Bessel term"""
    sqrt = np.sqrt(1 - 4*x**2)
    return (1 - 14*x**2 - 2*x**4 - 12*x**6)*sqrt + 12*x**4*(x**4 - 1)*L_func(x)

def f2(x):
    """Eq. (2.12), second Bessel term"""
    sqrt = np.sqrt(1 - 4*x**2)
    return 4*( x**2*(2 + 10*x**2 - 12*x**4)*sqrt + 6*x**4*(1 - 2*x**2 + 2*x**4)*L_func(x) )

def gamma_nu_ll(m_N, me, mixing_parameter_squared):
    """Eq. (2.10), N → ν_e ℓ⁻ ℓ⁺"""
    x = me/m_N
    C1 = 0.25*(1 - 4*sin2w + 8*sin2w**2)
    C2 = 0.5*(-sin2w + 2*sin2w**2)
    return mixing_parameter_squared * (fermi_constant**2 * m_N**5)/(96*np.pi**3) * ( (C1+2*sin2w)*f1(x) + (C2+sin2w)*f2(x) )

def gamma_e_l_nu(m_N, ml, mixing_parameter_squared):
    """Eq. (2.14), N → e⁻ ℓ⁺ νℓ"""
    x = ml/m_N
    return mixing_parameter_squared * (fermi_constant**2 * m_N**5)/(192*np.pi**3) * (1 - 8*x**2 + 8*x**6 - x**8 - 12*x**4*np.log(x**2))

def gamma_3nu(m_N, mixing_parameter_squared):
    """Eq. (2.15), N → ν ν ν"""
    return mixing_parameter_squared * (fermi_constant**2 * m_N**5)/(96*np.pi**3)

#Idea solve my running through the indexes of the array ie array[i]
def gamma_P_nu(m_N, mP, fP, mixing_parameter_squared):
    """Eq. (2.16), N → P⁰ ν"""
    xP = mP/m_N
    return mixing_parameter_squared * ((fermi_constant**2 * m_N**3 * fP**2)/(32*np.pi)) * (1 - xP**2)**2

def gamma_Pp_e(m_N, mP, fP, Vqq, mixing_parameter_squared):
    """Eq. (2.17), N → P⁺ e⁻"""
    xP = mP/m_N
    xe = 0  # assume m_e negligible
    lam = (1 - xP**2 - xe**2)**2 - 4*xP**2*xe**2
    return mixing_parameter_squared * (fermi_constant**2 * m_N**3 * fP**2 * Vqq**2)/(16*np.pi) * np.sqrt(lam) * ((1 - xP**2)**2 - xe**2*(2 + xP**2 - xe**2))

def gamma_V_nu(m_N, mV, fV, kappaV, mixing_parameter_squared):
    """Eq. (2.18), N → V⁰ ν"""
    xV = mV/m_N
    return mixing_parameter_squared * (fermi_constant**2 * m_N**3 * fV**2 * kappaV**2)/(32*np.pi * mV**2) * (1 + 2*xV**2)*(1 - xV**2)**2

def gamma_Vp_e(m_N, mV, fV, Vqq, mixing_parameter_squared):
    """Eq. (2.19), N → V⁺ e⁻"""
    xV = mV/m_N
    xe = 0
    lam = (1 - xV**2 - xe**2)**2 - 4*xV**2*xe**2
    return mixing_parameter_squared * (fermi_constant**2 * m_N**3 * fV**2 * Vqq**2)/(16*np.pi * mV**2) * np.sqrt(lam) * ((1 - xV**2)*(1 + 2*xV**2) + xe**2*(xV**2 + xe**2 - 2))

def gamma_N_to_SM(m_N, mixing_parameter_squared):
    """Eq. (2.20): sum over all open SM channels."""
    g = 0.
    # leptonic
    for mℓ in [0.000511] :  # lepton masses in GeV, due to the square root of the equation only the electron is possible 
        g += gamma_nu_ll(m_N, mℓ, mixing_parameter_squared)
    for mℓ in [0.1057, 1.777]:
        g += 2*gamma_e_l_nu(m_N, mℓ, mixing_parameter_squared)
    g += 3*gamma_3nu(m_N, mixing_parameter_squared) # for each of the 3 possible neutrino generations
    # pseudoscalars neutral, good trial
    pseudoscalar_neutral_mesons = [
        {'m':0.1349766, 'f':0.13},  # pion_neutral
        {'m':0.497614 , 'f':0.156},  # kaon neutral 
        {'m':0.547862 , 'f':0.1},   #eta  find decay constants eta, prime
        {'m':0.95778 , 'f':0.1},   #eta prime, this is an estimate
        
    ]
    for P in pseudoscalar_neutral_mesons:
        g += gamma_P_nu(m_N, P['m'], P['f'], mixing_parameter_squared)
        
        #g += 2*gamma_Pp_e(m_N, P['m'], P['f'], P['CKM'], mixing_parameter_squared)
    # vector mesons (example)
    pseudoscalar_charged_mesons = [
        {'m':0.13957, 'f':0.13, 'CKM':0.974}, #pion charged, CKM: |V_ud|
        {'m':0.493677, 'f':0.156, 'CKM': 0.225}, #kaon charged CKM: |V_us|
        {'m':1.86961, 'f':0.212, 'CKM':0.225}, # D meson charged |V_cd|
        {'m':1.96830, 'f':0.249, 'CKM':0.973}, # D strange meson charged 
    ]
    
    for P in pseudoscalar_charged_mesons:
        g += 2*gamma_Pp_e(m_N, P['m'], P['f'], P['CKM'], mixing_parameter_squared)
        # Couldn'y find a good source for the vector couplings so I just added test value
        ############################
        
    vector_neutral_mesons = [
        {'m':0.77526, 'f':0.171, 'kappa':0.038 }, #rho vector meson
        {'m':0.78265, 'f':0.155, 'kappa':0.038}, #omega vector meson
        {'m':1.019461, 'f':0.232, 'kappa': 0.038}, #phi vector meson
        {'m':0.89581, 'f':0.178,'kappa': 0.038}, # kaon_star_neutral meson
    ]
    
    for V in vector_neutral_mesons:
        g += gamma_V_nu(m_N, V['m'], V['f'], V['kappa'], mixing_parameter_squared)
        
    vector_charged_mesons = [
        {'m': 0.77511, 'f':0.171, 'CKM':0.974 },  #rho charged meson V_ud
        {'m':0.89166, 'f':0.178, 'CKM':0.225 },  #kaon_star_charged meson V_us

    ]
        
    for V in vector_charged_mesons:
        #g += gamma_V_nu(m_N, V['m'], V['f'], V['kappa'], mixing_parameter_squared)
        g += 2*gamma_Vp_e(m_N, V['m'], V['f'], V['CKM'], mixing_parameter_squared)
    return g



# quick sanity check at m_N=0.5 GeV:
gamma_N_SM = gamma_N_to_SM(0.1, mixing_parameter_squared)





def decay_width_N_av(f_a,m_N,mixing_parameter_squared, m_a):
    decay_width_N_av = mixing_parameter_squared*(m_N**3)*(1/(4*np.pi*(f_a**2)))*(np.sqrt(1+ (m_a/m_N)**2))*((1-(m_a/m_N)**2)**1.5)
    return decay_width_N_av

def decay_width_a_vv(f_a, m_N, mixing_parameter_squared, m_a, m_v):
    #decay_width_a_vv = (1/2* np.pi * f_a**2)*m_N**2 * m_a * mixing_parameter_squared**2 * (np.sqrt(1- (4*m_v**2)/m_a**2)) * (1- (2*m_v**2)/m_a**2)
    decay_width_a_vv = ((m_N**2 )* m_a * (mixing_parameter_squared**2))/(2* np.pi * f_a**2)
    return decay_width_a_vv
'''
# This is the part outside the main brackets
def cross_section_aa_NN(m_N,m_a,mixing_parameter_squared, s, f_a):
    
    # Check kinematic thresholds: s must be > max(4m_a², 4m_N²)
    #s_threshold = max(4 * m_a**2, 4 * m_N**2)
    s_threshold = 4.1 * m_a**2
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

    #acoth_val = 0.5 * np.log((acoth_arg + 1) / (acoth_arg - 1))
    acoth_val = mpmath.acoth(acoth_arg)
    term_C = (2 * numerator_C * acoth_val) / denominator_C

    # Combine terms B and C
    bracket_term = term_B + term_C

    # Prefactor
    prefactor = (4 * m_N**4 * sqrt_s_minus_4ma2 ) / (np.pi * f_a**4 * s**1.5)

    # Full cross-section
    cross_section = prefactor * bracket_term
    return cross_section
'''
'''
def cross_section_aa_NN(m_N,m_a,mixing_parameter_squared, s, f_a):
    
    # Check kinematic thresholds
    s_threshold = 4.1 * m_a**2
    if s <= s_threshold:
        return 0.0
    
    sqrt_s_minus_4ma2 = np.sqrt(s - 4 * m_a**2)
    sqrt_s_minus_4mN2 = np.sqrt(s - 4 * m_N**2)
    prefactor_numerator = 4*m_N**2 *sqrt_s_minus_4ma2
    prefactor_denominator = np.pi*f_a**4 *s**1.5
    prefactor = prefactor_numerator/prefactor_denominator
    termA_numerator = m_a**4 - 4*m_a**2 *m_N**2
    termA_denominator = m_a**4 - 4*m_a**2 *m_N**2+ m_N**2 *s
    termA = termA_numerator/termA_denominator
    acoth_term_numerator  =  s-2*m_a**2
    acoth_term_denominator = sqrt_s_minus_4ma2*sqrt_s_minus_4mN2
    acoth_term = mpmath.acoth(acoth_term_numerator/acoth_term_denominator)
    termB_numerator = 2*(s-2*m_a**2)*acoth_term 
    termB_denominator = sqrt_s_minus_4ma2*sqrt_s_minus_4mN2
    termB = termB_numerator/termB_denominator
    cross_section = prefactor*(-1-termA+termB)
    return cross_section
'''
def cross_section_aa_NN(m_N,m_a,mixing_parameter_squared, s, f_a):
    
    # Check kinematic thresholds
    s_threshold = 4.1 * m_a**2
    if s <= s_threshold:
        return 0.0
    
    sqrt_s_minus_4ma2 = np.sqrt(s - 4 * m_a**2)
    sqrt_s_minus_4mN2 = np.sqrt(s - 4 * m_N**2)
    
    prefactor_numerator = 4*m_N**4 *sqrt_s_minus_4ma2 * mixing_parameter_squared**2
    prefactor_denominator = np.pi*f_a**4 *s**1.5
    prefactor = prefactor_numerator/prefactor_denominator
    termA_numerator = -(2*m_a**4 - 4*m_a**2 *m_N**2 + 4*m_N**4 -m_N**2 *s)
    termA_denominator = m_a**4 - 2*m_a**2 * m_N**2 + m_N**4
    termA = termA_numerator/termA_denominator
    acoth_numerator = 2*(m_a**2 + m_N**2) - s 
    acoth_denominator = sqrt_s_minus_4ma2*sqrt_s_minus_4mN2
    acoth_term = mpmath.acoth(acoth_numerator/acoth_denominator)
    termB_numerator = 2*(2*(m_a**2 + m_N**2)- s)*acoth_term 
    termB_denominator = sqrt_s_minus_4ma2*sqrt_s_minus_4mN2
    termB = termB_numerator/termB_denominator 
    cross_section = prefactor*(termA + termB)
    return cross_section

    

def cross_section_aa_vN(m_N, m_a,mixing_parameter_squared,s, f_a): #A.7 of paper
    
    s_threshold = 4.1 * m_a**2
    if s <= s_threshold:
        return 0.0

    sqrt_s_minus_4ma2 = np.sqrt(s - 4 * m_a**2)
    sqrt_s = np.sqrt(s)
    # prefactor
    prefactor_numerator = 2*m_N**4 * sqrt_s_minus_4ma2 * mixing_parameter_squared**2
    prefactor_denominator = np.pi*f_a**4 *s**1.5
    prefactor = prefactor_numerator/prefactor_denominator
    termA_numerator = -2*(2*m_a**4 - 7*m_a**2 * m_N**2 + 5*m_N**4 + m_N**2 *s)
    termA_denominator = m_a**4 - 2*m_a**2 *m_N**2 + m_N**4 + m_N**2 *s
    termA = termA_numerator/ termA_denominator
    coth_term_numerator = 2*(m_a**2 - m_N**2) - s
    coth_term_denominator = sqrt_s_minus_4ma2 *sqrt_s
    coth_term = mpmath.acoth(coth_term_numerator/coth_term_denominator)
    termB_numerator = 4*(2*m_a**2 - 4*m_N**2 -s)*coth_term
    termB_denominator = sqrt_s_minus_4ma2 * sqrt_s
    termB = termB_numerator/termB_denominator
    cross_section = prefactor*(termA + termB)
    return cross_section



def cross_section_aa_vN_v_med(m_N, m_a,mixing_parameter_squared,s, f_a): #A.10 of paper
    
    s_threshold = 4.1 * m_a**2
    if s <= s_threshold:
        return 0.0
    
    sqrt_s_minus_4ma2 = np.sqrt(s - 4 * m_a**2)
    sqrt_s = np.sqrt(s)
    prefactor_numerator = 4*m_N**4*sqrt_s_minus_4ma2*mixing_parameter_squared**3 
    prefactor_denominator = np.pi* f_a**4 *s**1.5
    prefactor = prefactor_numerator / prefactor_denominator 
    termA = (m_N**2 / m_a**2)
    acoth_term_numerator = s- 2*m_a**2
    acoth_term_denominator = sqrt_s_minus_4ma2* sqrt_s
    acoth_term = mpmath.acoth(acoth_term_numerator/ acoth_term_denominator)
    termB_numerator = 2*(s- 2*m_a**2 )*acoth_term 
    termB_denominator = sqrt_s_minus_4ma2 *sqrt_s 
    termB = termB_numerator / termB_denominator 
    cross_section = prefactor * (-2 + termA + termB)
    
    return cross_section 

    
    
    
    


from scipy.integrate import quad
from scipy.special import kn, k1

def gamma_aa_NN(z, m_N, m_a, mixing_parameter_squared, f_a):

   
    #s_min = max((4.1*m_a)**2, (4.1*m_N)**2)
    s_min = s_threshold = 4 * m_a**2
    
    
    s_max = np.inf
    

    def lambda_func(a, b, c):
        return (a - b - c)**2 - 4*b*c
    

    def reduced_sigma(s):

        sigma_val = cross_section_aa_NN(m_N, m_a, mixing_parameter_squared, s, f_a)
        
        #print(sigma_val, 'sigma_val')
        lam = lambda_func(1, m_a**2/s, m_a**2/s)
        #print(lam, 'lam')
  
        return 2 * s * sigma_val * lam
    

    def integrand(s):
        # Handle snear zero to avoid numerical issues
        if s <= 0:
            return 0.0
        return np.sqrt(s) * reduced_sigma(s) * k1(np.sqrt(s)/(m_N/z))
    

    # below format ignores error value as integrand produces result and error
    integral, _ = quad(integrand, s_min, s_max)
    print(integral, 'integral') # problem with this integral

    gamma = ((m_N/z) / (64 * np.pi**4)) * integral
    return gamma


def gamma_aa_vN(z, m_N, m_a, mixing_parameter_squared, f_a):

   
    #s_min = max((4.1*m_a)**2, (4.1*m_N)**2)
    s_min = 4.1*m_a**2

    #s_max = (100*(m_N/z))**2  # Safe upper bound where K₁ decays sufficiently
    s_max = np.inf
    

    def lambda_func(a, b, c):
        return (a - b - c)**2 - 4*b*c
    

    def reduced_sigma(s):

        sigma_val = cross_section_aa_vN(m_N, m_a, mixing_parameter_squared, s, f_a)
        

        lam = lambda_func(1, m_a**2/s, m_a**2/s)
  
        return 2 * s * sigma_val * lam
    

    def integrand(s):
        # Handle snear zero to avoid numerical issues
        if s <= 0:
            return 0.0
        return np.sqrt(s) * reduced_sigma(s) * k1(np.sqrt(s)/(m_N/z))
    

    # below format ignores error value as integrand produces result and error
    integral, _ = quad(integrand, s_min, s_max, 
                      epsabs=1e-9, epsrel=1e-6, 
                      limit=1000)
    

    gamma = ((m_N/z) / (64 * np.pi**4)) * integral
    return gamma


'''
#A.10
def gamma_aa_vN_v_med(z, m_N, m_a, mixing_parameter_squared, f_a):

   
    #s_min = max((4.1*m_a)**2, (4.1*m_N)**2)
    s_min = 4.1*m_a**2

    #s_max = (100*(m_N/z))**2  # Safe upper bound where K₁ decays sufficiently
    s_max = np.inf
    

    def lambda_func(a, b, c):
        return (a - b - c)**2 - 4*b*c
    

    def reduced_sigma(s):

        sigma_val = cross_section_aa_vN_v_med(m_N, m_a, mixing_parameter_squared, s, f_a)
        

        lam = lambda_func(1, m_a**2/s, m_a**2/s)
  
        return 2 * s * sigma_val * lam
    

    def integrand(s):
        # Handle snear zero to avoid numerical issues
        if s <= 0:
            return 0.0
        return np.sqrt(s) * reduced_sigma(s) * k1(np.sqrt(s)/(m_N/z))
    

    # below format ignores error value as integrand produces result and error
    integral, _ = quad(integrand, s_min, s_max, 
                      epsabs=1e-9, epsrel=1e-6, 
                      limit=1000)
    

    gamma = ((m_N/z) / (64 * np.pi**4)) * integral
    return gamma

'''


#implementing 3.3 of new paper, where K refer to the Bessel functions.
def thermally_averaged_decay_width_density_N_SM(n_N_eq, decay_width_N_SM,z):
    decay_width_N_SM = gamma_N_to_SM(0.1, mixing_parameter_squared)
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
gamma_N_SM_n_N_eq_ratio_list =  []
gamma_N_av_n_N_eq_ratio_list = []
gamma_a_vv_n_N_eq_ratio_list = []

gamma_aa_NN_n_a_eq_list = []

gamma_aa_vN_n_a_eq_list = []
gamma_aa_vN_v_med_n_a_eq_list = []


for z in mass_temp_ratio:
    if z > 0.1/173.3 and z < : # !!! define inequality for all of these!!!
        g_star = 106.75
    elif z > 0.1/125.6:
        g_star = 96.25
    elif z> 0.1/91.2:
        g_star = 95.25
    elif z> 0.1/80.4:
        g_star = 92.25
    elif z> 0.1/4.19:
        g_star = 86.25
    elif z>0.1/1.777:
        g_star = 75.75
    elif z> 0.1/1.29:
        g_star = 72.25
    elif z> 0.1/0.214:
        g_star = 61.75
    elif z> 0.1/0.1396:
        g_star = 17.25
    elif z> 0.1/0.135:
        g_star = 15.25
    elif z> 0.1/0.1057:
        g_star = 14.25
    elif z>0.1/800e-6:
        g_star = 10.75
    #for the last two used g*s, come back and check this 
    elif z>0.1/511e-3:
        g_star = 0.1/7.409
    else:
        g_star = 3.909
    
    
    H_value = H(z)
    s_value = entropy_density(z)
    n_N_eq_value = equilibrium_species_density(2, 0.1, z)
    Y_N_eq_value = n_N_eq_value / entropy_density(z)
    #Now ALPs
    n_a_eq_value = equilibrium_species_density(2, 1e-6,z)
    decay_width_N_SM = gamma_N_to_SM(0.1, mixing_parameter_squared)
    gamma_value = thermally_averaged_decay_width_density_N_SM(n_N_eq_value, decay_width_N_SM, z)
    #Below I am calculating the ratio of 
    gamma_N_SM_n_N_eq_ratio_value = decay_width_N_SM *(kn(1,z)/kn(2,z))
    gamma_N_av_n_N_eq_ratio_value =  decay_width_N_av(1e3,m_N,mixing_parameter_squared, m_a)*(kn(1,z)/kn(2,z))
    gamma_a_vv_n_N_eq_ratio_value = decay_width_a_vv(f_a, m_N, mixing_parameter_squared, m_a, m_v)*(kn(1,z)/kn(2,z))
    
    gamma_aa_NN_value = gamma_aa_NN(z, m_N, m_a, mixing_parameter_squared, f_a)
    gamma_aa_NN_n_a_eq_ratio_value = gamma_aa_NN_value/n_a_eq_value
    
    gamma_aa_vN_value = gamma_aa_vN(z, m_N, m_a, mixing_parameter_squared, f_a)
    gamma_aa_vN_n_a_eq_ratio_value = gamma_aa_vN_value/n_a_eq_value
    #A,10
    #gamma_aa_vN_v_med_value = gamma_aa_vN_v_med(z, m_N, m_a, mixing_parameter_squared, f_a)
    #gamma_aa_vN_v_med_n_a_eq_ratio_value = gamma_aa_vN_v_med_value/n_a_eq_value
    
    

     #adding to the lists of values
    H_list.append(H_value)
    s_list.append(s_value)
    n_N_eq_list.append(n_N_eq_value)
    Y_N_eq_list.append(Y_N_eq_value)
    gamma_list.append(gamma_value)
    gamma_N_SM_n_N_eq_ratio_list.append(gamma_N_SM_n_N_eq_ratio_value)
    gamma_N_av_n_N_eq_ratio_list.append(gamma_N_av_n_N_eq_ratio_value)
    gamma_aa_NN_n_a_eq_list.append(gamma_aa_NN_n_a_eq_ratio_value)
    gamma_a_vv_n_N_eq_ratio_list.append(gamma_a_vv_n_N_eq_ratio_value)
    gamma_aa_vN_n_a_eq_list.append(gamma_aa_vN_n_a_eq_ratio_value)
   # gamma_aa_vN_v_med_n_a_eq_list.append(gamma_aa_vN_v_med_n_a_eq_ratio_value) A.10
    
#was here

plt.loglog(mass_temp_ratio, gamma_N_av_n_N_eq_ratio_list, label = r'$\gamma_{N \rightarrow a \nu }/ n_{N}^{eq}$')
plt.loglog(mass_temp_ratio, gamma_N_SM_n_N_eq_ratio_list, label = r'$\gamma_{N \rightarrow SM }/ n_{N}^{eq}$')
plt.loglog(mass_temp_ratio, H_list , label = 'Hubble Parameter')
plt.loglog(mass_temp_ratio, gamma_a_vv_n_N_eq_ratio_list,label = r'$\gamma_{a \rightarrow \nu \nu }/ n_{a}^{eq}$' )
#plt.loglog(mass_temp_ratio, gamma_aa_NN_n_a_eq_list)
plt.loglog(mass_temp_ratio,gamma_aa_NN_n_a_eq_list, label = 'aa->NN')
plt.loglog(mass_temp_ratio, gamma_aa_vN_n_a_eq_list, label = 'aa ->vN (N)')
#plt.loglog(mass_temp_ratio, gamma_aa_vN_v_med_n_a_eq_list, label = 'aa ->vN (v)') A.10
plt.xlim(1e-4, 1e2)
plt.ylim(1e-45, max(H_list))
plt.title(r'$\gamma_{X \rightarrow Y }/ n_{X}^{eq}$')
plt.xlabel("z = m_N /T")
plt.ylabel(r'$\gamma_{X \rightarrow Y }/ n_{X}^{eq}$')
plt.legend()
plt.show()
print("End")
#need to fix the definitions of z

    
