#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import kn, k1, zeta
import matplotlib.pyplot as plt

# =============================================================================
# 1. Physical Constants
# =============================================================================
# Using natural units (GeV)
M_pl = 2.435e18      # Reduced Planck mass [GeV]
fermi_constant = 1.166e-5 # Fermi constant [GeV^-2]
g_star = 10.75        # Relativistic degrees of freedom for the relevent pre-BBN period

# Model Parameters
m_N = 0.1             # HNL Mass [GeV]
m_a = 1e-6            # ALP Mass [GeV]
m_nu = 1e-9           # Active neutrino mass [GeV]
mixing_parameter = 1e-10 # Active-sterile neutrino sector mixing
f_a = 1e12            # ALP decay constant [GeV]
f_M = 0.1             # Meson decay constant [GeV]
# Integration limits
z_init = 0.3675
z_final = 120

# =============================================================================
# 2. Cosmological Thermodynamic Functions
# =============================================================================
def H(z, m_N_val, M_pl_val):
    """Calculates the Hubble rate H(z)."""
    T = m_N_val / z
    return np.sqrt((np.pi**2 * g_star) / 90) * (T**2 / M_pl_val)

def entropy_density(z, m_N_val, g_star_val):
    """Calculates the entropy density s(z)."""
    T = m_N_val / z
    return ((2 * (np.pi**2)) / 45) * g_star_val * (T**3)

def equilibrium_yield(z, m, g, p):
    """
    Calculates the equilibrium comoving yield Y_eq = n_eq / s.
    Handles both relativistic (m << T) and non-relativistic (m >> T) cases.
    """
    T = p['m_N'] / z
    x = m / T

    if np.any(x < 0.1):  # Particle is relativistic
        if g == 1: # Boson
            n_eq = (zeta(3) / (np.pi**2)) * g * T**3
        else: # Fermion
            n_eq = (3/4) * (zeta(3) / (np.pi**2)) * g * T**3
    else:  # Particle is non-relativistic
        n_eq = g * (m * T / (2 * np.pi))**1.5 * np.exp(-x)

    s_z = entropy_density(z, p['m_N'], p['g_star'])
    return n_eq / s_z if s_z > 0 else 0

def equilibrium_species_density(g,mass,z):
    species_density = g*(mass**3)*(((1/(2*np.pi*z))**1.5))*np.exp(-z)
    return species_density 
n_N_eq = equilibrium_species_density(2,0.1,z_init)    
n_a_eq = equilibrium_species_density(1,m_a,z_init)


# =============================================================================
# 3. Interaction Rates
# =============================================================================
def decay_width_N_SM(f_M_val, m_N_val, mix_val):
    term1 = (6 * f_M_val**2) / np.pi
    term2 = m_N_val**2 / (20 * np.pi**3)
    return (term1 + term2) * (m_N_val**3) * (abs(mix_val)**2) * (fermi_constant**2)

def decay_width_N_av(mix_val, m_N_val, f_a_val):
    return ((abs(mix_val)**2) * (m_N_val**3)) / (4 *np.pi*f_a_val**2)

def decay_width_a_vv(mix_val, m_N_val, f_a_val, m_a_val):
    return ((m_N_val**2) * m_a_val * (mix_val**4)) / (2 * np.pi * (f_a_val**2))

# =============================================================================
# 4. The Coupled Boltzmann Equations (DEBUGGING VERSION)
# =============================================================================
def coupled_boltzmann_eqs(z, Y, p):

    try:
        # This print statement will show us every time the solver calls this function
        print(f"--- ODE function called with z = {z}")

        Y_N, Y_a = Y
        T = p['m_N'] / z

        Y_N_eq = equilibrium_yield(z, p['m_N'], 2, p)
        Y_a_eq = equilibrium_yield(z, p['m_a'], 1, p)


        def bessel_ratio(m, T_val):
            z_species = m / T_val
            # Use np.any/np.divide for safety in case solver uses vector z
            #below if statement is for optimization, to avoid the need for calculations
            # for bessel ratio k1/k2 very close to one. Also avoids issues with precision values
            
            if np.any(z_species < 1e-3):
                if isinstance(z_species, np.ndarray):
                    return np.ones_like(z_species)
                return 1.0
            k1_val = kn(1, z_species)
            k2_val = kn(2, z_species)
            return np.divide(k1_val, k2_val, out=np.zeros_like(k1_val, dtype=float), where=k2_val!=0)

        Gamma_N_SM_avg = p['Gamma_N_SM'] * bessel_ratio(p['m_N'], T)*n_N_eq
        Gamma_N_av_avg = p['Gamma_N_av'] * bessel_ratio(p['m_N'], T)
        Gamma_a_vv_avg = p['Gamma_a_vv'] * bessel_ratio(p['m_a'], T)
        print("    ...done.")

        # Annihilation terms are off
        sv_aa_NN = 0.0
        sv_aa_vv = 0.0

        print("    Calculating s(z) and H(z)...")
        s_z = entropy_density(z, p['m_N'], p['g_star'])
        H_z = H(z, p['m_N'], p['M_pl'])
        prefactor = 1 / (z * H_z * s_z)
        print("    ...done.")

        print("    Calculating collision terms...")
        term_N1 = 0.0 # Annihilation is off
        term_N2 = s_z * Gamma_N_SM_avg * (Y_N - Y_N_eq)
        term_N3 = s_z * Gamma_N_av_avg * (Y_N - (Y_N_eq / Y_a_eq) * Y_a) if Y_a_eq > 0 else 0
        dY_N_dz = -prefactor * (term_N1 + term_N2 + term_N3)

        term_a1 = 0.0 # Annihilation is off
        term_a2 = -term_N1
        term_a3 = s_z * Gamma_a_vv_avg * (Y_a - Y_a_eq)
        term_a4 = -term_N3
        dY_a_dz = -prefactor * (term_a1 + term_a2 + term_a3 + term_a4)
        print("    ...done. Returning derivatives.")

        return [dY_N_dz, dY_a_dz]

    # Was having issues with errors in the differential equations at low z values
     # Attempted to solve this by changing the integration method to 'Radau' which
     # is better suited to stiffer equations, this if and except statement was used
     # to test a number of methods 
    except Exception as e:
        print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!!!   ERROR in coupled_boltzmann_eqs    !!!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(f"The error happened when the solver was at z = {z}")
        print(f"The error message is: {e}")
        print(f"The type of the error is: {type(e)}")
        raise # Re-raise the error to stop the program as before


# =============================================================================
# 5. Solve the Equations
# =============================================================================
# Pre-calculate constant decay widths
Gamma_N_SM_val = decay_width_N_SM(f_M, m_N, mixing_parameter)
Gamma_N_av_val = decay_width_N_av(mixing_parameter, m_N, f_a)
Gamma_a_vv_val = decay_width_a_vv(mixing_parameter, m_N, f_a, m_a)

# Package all parameters into a dictionary for clean passing
params = {
    'm_N': m_N, 'm_a': m_a, 'mix': mixing_parameter, 'f_a': f_a, 'f_M': f_M,
    'M_pl': M_pl, 'g_star': g_star,
    'Gamma_N_SM': Gamma_N_SM_val, 'Gamma_N_av': Gamma_N_av_val, 'Gamma_a_vv': Gamma_a_vv_val
}



# Initial conditions
Y_N_initial = 1e-40
Y_a_initial = 1e-40
Y0 = [Y_N_initial, Y_a_initial]


# Solve the coupled ODEs
sol = solve_ivp(
    fun=lambda z, Y: coupled_boltzmann_eqs(z, Y, params),
    t_span=[z_init, z_final],
    y0=Y0,
    method='Radau',
    dense_output=True,
    atol=1e-40,
    rtol=1e-7
)

# This part will only be reached if the solver SUCCEEDS
print("No Crash")
print(sol)


# =============================================================================
# 6. Plot the Results
# =============================================================================

# Create a log-spaced array for plotting
z_plot = np.logspace(np.log10(z_init), np.log10(sol.t[-1]), 500)

# Evaluate the solution at the plotting points
Y_sol = sol.sol(z_plot)
Y_N_vals = Y_sol[0]
Y_a_vals = Y_sol[1]

# Also calculate the equilibrium yields for comparison
Y_N_eq_vals = np.array([equilibrium_yield(z, m_N, 2, params) for z in z_plot])
Y_a_eq_vals = np.array([equilibrium_yield(z, m_a, 1, params) for z in z_plot])

# --- Plotting ---
plt.figure(figsize=(12, 8))

plt.loglog(z_plot, Y_N_vals, label=r'$Y_N$ (HNL)', color='orange', linewidth=2.5)
plt.loglog(z_plot, Y_a_vals, label=r'$Y_a$ (ALP)', color='blue', linewidth=2.5)
plt.loglog(z_plot, Y_N_eq_vals, label=r'$Y_N^{\mathrm{eq}}$', color='orange', linestyle='--', alpha=0.7)
plt.loglog(z_plot, Y_a_eq_vals, label=r'$Y_a^{\mathrm{eq}}$', color='blue', linestyle='--', alpha=0.7)

plt.xlabel(r'Inverse Temperature $z = m_N / T$', fontsize=14)
plt.ylabel('Comoving Yield $Y = n/s$', fontsize=14)
plt.title('Evolution of HNL and ALP Abundances', fontsize=16)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend(fontsize=12)

# Dynamically set y-axis limits to focus on the solved curves
y_min_val = 1e-30 # Set a floor for the y-axis
y_max_val = max(np.max(Y_N_vals), np.max(Y_a_vals), np.max(Y_a_eq_vals))
plt.ylim(y_min_val, y_max_val * 10)
plt.xlim(z_init, z_final)

plt.show()

print("\n--- Script finished. ---")