#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 15:13:19 2025

@author: janwawrzynek
"""
# consider 10e-10 vs 1e-10
# attemoting 3.8 and 3.9
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import kn
import matplotlib.pyplot as plt
import mpmath

# Physical constants (natural units: GeV, s^-1)
M_pl = 2.435e18          # Planck mass [GeV]
fermi_constant = 1.166e-5           # Fermi constant [GeV^-2]
m_N = 1e-1 # HNL Mass in GeV
m_a = 1e-6
m_nu = 1e-9 #mass 1eV for now
mixing_parameter = 1e-10 # active-sterile neutrino sector mixing
f_a = 1e12 # ALP decay constant scenario 1
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
n_a_eq = equilibrium_species_density(1,m_a,z_init)
    
#def comoving_yield(species_density,entropy_density):
#    comoving_yield = species_density/entropy_density
#    return comoving_yield
#HNL comoving yield at thermal equilibrium Y_N_eq
#Y_N_eq = comoving_yield(n_N_eq,entropy_density)
# above didn't work
Y_N_eq = n_N_eq / entropy_density(z_init)
Y_a_eq = n_a_eq/entropy_density(z_init)

    
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
gamma_N_SM = thermally_averaged_decay_width_density_N_SM(n_N_eq, decay_width_N_SM,z)


#Defining the dominant decay width of N_> av
def decay_width_N_av(mixing_parameter, M_N, f_a):
    decay_width_N_av = ((abs(mixing_parameter)**2)*(M_N**3))/4*np.pi**2
    return decay_width_N_av
#f_a = 1 TeV =  1e+3 GeV
decay_width_N_av =  decay_width_N_av(1-10,0.1,1e+3)
def thermally_averaged_decay_width_density_N_av(n_a_eq, decay_width_N_av,z):
    thermally_averaged_interaction_density_N_av = n_a_eq*(kn(1,z_vals)/kn(2,z_vals))*decay_width_N_av
    return thermally_averaged_interaction_density_N_av
gamma_N_av = thermally_averaged_decay_width_density_N_av(n_a_eq, decay_width_N_av,z)



def decay_width_a_vv(mixing_parameter, m_N,f_a,m_a):
    decay_width_a_vv = ((m_N**2)*m_a*(mixing_parameter**4))/(2*np.pi*(f_a**2))
    return decay_width_a_vv
decay_width_a_vv = decay_width_a_vv(1e-10, 0.1,1e+3,1e-6)

def thermally_averaged_decay_width_density_a_vv(n_a_eq, decay_width_a_vv,z):
    thermally_averaged_interaction_density_a_vv = n_a_eq*(kn(1,z_vals)/kn(2,z_vals))*decay_width_a_vv
    return thermally_averaged_interaction_density_a_vv
gamma_a_vv = thermally_averaged_decay_width_density_a_vv(n_a_eq, decay_width_a_vv,z)



# Calculating the cross section of aa->NN using equation  A.8 
#for now consider the minimum mandelstam s variable value s> 2m_N**2
def mandelstam_s(mass):
    mandelstam_s = 2*(mass**2)
    return mandelstam_s

s_value = mandelstam_s(0.1)

    


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

    # Term C: acoth term
    numerator_C = s - 2 * m_a**2
    denominator_C = sqrt_s_minus_4ma2 * sqrt_s_minus_4mN2
    acoth_arg = numerator_C / denominator_C
    # Check domain for acoth (|arg| must be > 1)
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



s_physical =  4.1*m_N**2
cross_section_val_aa_NN = cross_section_aa_NN(0.1,1e-6,1e-10, s_physical, 10e12)

def cross_section_aa_vv(m_N, m_a, mixing_parameter, s, f_a):
    # Check kinematic threshold: s must be > 4m_a²
    if s <= 4 * m_a**2:
        return 0.0

    # Precompute common terms
    sqrt_s = np.sqrt(s)
    s_minus_4ma2 = s - 4 * m_a**2
    sqrt_s_minus_4ma2 = np.sqrt(s_minus_4ma2)
    
    # Argument for coth⁻¹ function
    numerator_arg = s - 2 * m_a**2
    denominator_arg = sqrt_s_minus_4ma2 * sqrt_s
    arg = numerator_arg / denominator_arg

    # Check domain for acoth (|arg| must be > 1)
    if abs(arg) <= 1:
        return 0.0

    # Calculate coth⁻¹(x) = 0.5 * ln((x+1)/(x-1))
    acoth_val = 0.5 * np.log((arg + 1) / (arg - 1))

    # Compute the term in brackets
    bracket_term = -2 + (2 * (s - 2 * m_a**2) * acoth_val) / (sqrt_s_minus_4ma2 * sqrt_s)
    # Full cross-section calculation
    cross_section = (
        (4 * m_N**4 * sqrt_s_minus_4ma2 * mixing_parameter**8) /
        (np.pi * f_a**4 * s**1.5)
    ) * bracket_term

    return cross_section
cross_section_val_aa_vv = cross_section_aa_vv(0.1,1e-6,1e-10, s_physical, 10e12)


from scipy.integrate import quad
from scipy.special import kn, k1

def gamma_aa_NN(T, m_N, m_a, mixing_parameter, f_a):
    """
    Calculate the thermally averaged cross-section density γ for aa <-> NN process
    using Equation 3.4.
    
    Parameters:
    T: float             # Temperature [GeV]
    m_N: float           # Heavy neutrino mass [GeV]
    m_a: float           # ALP mass [GeV]
    mixing_parameter: float  # Mixing parameter
    f_a: float           # ALP decay constant [GeV]
    
    Returns:
    gamma: float         # Thermally averaged cross-section density [GeV⁴]
    """

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
        return np.sqrt(s) * reduced_sigma(s) * k1(np.sqrt(s)/(T))
    

    # below format ignores error value as integrand produces result and error
    integral, _ = quad(integrand, s_min, s_max, 
                      epsabs=1e-10, epsrel=1e-6, 
                      limit=1000)
    

    gamma = (T) / (64 * np.pi**4) * integral
    return gamma

gamma_vals_aa_NN_list = []

# Loop over each temperature in the T array
for temp in T:
    # Calculate gamma for a single temperature
    gamma_val = gamma_aa_NN(temp, m_N, m_a, mixing_parameter, f_a)
    gamma_vals_aa_NN_list.append(gamma_val)

# Convert the list of results into a NumPy array
gamma_vals_aa_NN = np.array(gamma_vals_aa_NN_list)


def gamma_aa_vv(T, m_N, m_a, mixing_parameter, f_a, m_nu):
    
    # 1. Calculate s_min = Max[(2m_a)², (2m_N)²]
    s_min = max((2*m_a)**2, (2*m_nu)**2)
    
    # 2. Define upper limit for integration (using exponential decay of K₁)
    s_max = (10*(T))**2  # Safe upper bound where K₁ decays sufficiently
    
   
    def lambda_func(a, b, c):
        return (a - b - c)**2 - 4*b*c
    

    def reduced_sigma(s):
        # Calculate the cross-section σ(s)
        sigma_val = cross_section_aa_vv(m_N, m_a, mixing_parameter, s, f_a)
        
       
        lam = lambda_func(1, m_a**2/s, m_a**2/s)
        
        
        return 2 * s * sigma_val * lam
    

    def integrand(s):
        # Handle s near zero to avoid numerical issues
        if s <= 0:
            return 0.0
        return np.sqrt(s) * reduced_sigma(s) * k1(np.sqrt(s)/(T))
    
    
    # below format ignores error value as integrand produces result and error
    integral, _ = quad(integrand, s_min, s_max, 
                      epsabs=1e-10, epsrel=1e-6, 
                      limit=1000)
    

    gamma = (T) / (64 * np.pi**4) * integral
    return gamma



gamma_vals_aa_vv_list = []


# Loop over each temperature in the T array
for temp in T:
    # Calculating gamma for a single temperature
    gamma_val = gamma_aa_vv(temp, m_N, m_a, mixing_parameter, f_a, m_nu)
    gamma_vals_aa_vv_list.append(gamma_val)



gamma_vals_aa_vv = np.array(gamma_vals_aa_vv_list)




# 6. The Coupled Boltzmann Equations 3.8 and 3.9
# =============================================================================
def coupled_boltzmann_eqs(z_vals, Y, H_value, entropy_density_value):
    """
    Returns the derivatives [dY_N/dz, dY_a/dz] for the ODE solver.
    Y[0] = Y_N (HNL yield)
    Y[1] = Y_a (ALP yield)
    """
    Y_N, Y_a = Y[0], Y[1]
    T = m_N / z

    # --- Pre-calculate equilibrium yields at current z ---
    Y_N_eq_z = Y_N_eq(z, m_N, 2)
    Y_a_eq_z = Y_a_eq(z, m_a, 1)

    # --- Calculate all gamma terms at current temperature T ---



    #interpolated values to allow for solving the ode using solve_ivp
    H_interp = interp1d(z_vals, H_vals, kind='cubic', fill_value="extrapolate")
    s_interp = interp1d(z_vals, s_vals, kind='cubic', fill_value="extrapolate")
    n_N_eq_interp = interp1d(z_vals, n_N_eq_vals, kind='cubic', fill_value="extrapolate")
    entropy_density_value_interp = interp1d(z_vals, entropy_density_value, kind='cubic', fill_value="extrapolate")
    
    
    
    # Prefactor for the ODEs
    prefactor = -1 / (z * H_interp * entropy_density_value_interp)
    # HNL coupled equation first
    #term_N1 = come back to this once get answer
    term_N2 = gamma_vals_aa_NN*((Y_a/Y_a_eq)**2 - (Y_N/Y_N_eq)**2)
    term_N3 = -gamma_N_SM*((Y_N/Y_N_eq) -1 )
    term_N4 = -gamma_N_av*((Y_N/Y_N_eq)- (Y_a/Y_a_eq))
    
    dY_N_dz = prefactor*(term_N2 + term_N3 + term_N4)
    
    term_a1 = -gamma_aa_vv*((Y_a/Y_a_eq)**2 -1)
    term_a2 = -gamma_vals_aa_vv*((Y_a/Y_a_eq)**2 - (Y_N/Y_N_eq)**2 )
    term_a3 = -gamma_a_vv*((Y_a/Y_a_eq) -1)
    term_a4 = gamma_N_av*((Y_N/Y_N_eq)- (Y_a/Y_a_eq))
    
    

    dY_a_dz = prefactor*(term_a1+term_a2+term_a3+term_a4)
    
    return [dY_N_dz, dY_a_dz]

Y_N_initial = Y_N_eq
Y_a_initial = 1e-40



# Initial conditions
Y0 = [Y_N_initial, Y_a_initial] 

# Correct call to the solver
sol = solve_ivp(
    fun=coupled_boltzmann_eqs,
    t_span=[z_init, z_final],       
    y0=Y0,                         
    method='Radau',                 
    t_eval=z_vals,                                  
)




# =============================================================================
# 8. Plotting
# =============================================================================
z_plot = np.logspace(np.log10(z_init), np.log10(z_final), 500)
Y_sol = sol.sol(z_plot)
Y_N_vals = Y_sol[0]
Y_a_vals = Y_sol[1]

# Also plot the equilibrium yields for comparison
Y_N_eq_vals = Y_N_eq(z_plot, m_N, 2)
Y_a_eq_vals = Y_a_eq(z_plot, m_a, 1)

plt.figure(figsize=(10, 7))
plt.loglog(z_plot, Y_N_vals, label=r'$Y_N$ (HNL)', color='orange')
plt.loglog(z_plot, Y_a_vals, label=r'$Y_a$ (ALP)', color='blue')
plt.loglog(z_plot, Y_N_eq_vals, label=r'$Y_N^{eq}$', color='orange', linestyle='--')
plt.loglog(z_plot, Y_a_eq_vals, label=r'$Y_a^{eq}$', color='blue', linestyle='--')

plt.xlabel(r'$z = m_N / T$')
plt.ylabel('Comoving Yield Y')
plt.title('Evolution of HNL and ALP Abundances')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()
plt.ylim(1e-30, 1e2)
plt.xlim(z_init, z_final)
plt.show()








# Solving for Y_N using differential equation 
 

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
    H_z = float(H_interp(z))
    s_z = float(s_interp(z))
    n_N_eq_z = float(n_N_eq_interp(z))
    gamma_z = float(gamma_interp(z))
    return -(1 / (z * H_z * s_z)) * gamma_z * ((Y_N / n_N_eq_z) - 1)

# now solve boltzmann equation using scipy ivp 

sol = solve_ivp(dY_N_dz, [z_init, z_final], [Y_N_eq], dense_output=True)
z_vals = sol.t
Y_vals = sol.y[0]

plt.plot(z_vals,Y_vals)
plt.xlim(min(z_vals), max(z_vals))
plt.xlabel("z = m_N /T")
plt.show()

#need to fix the definitions of z
