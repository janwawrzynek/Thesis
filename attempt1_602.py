#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 23:03:41 2025

@author: janwawrzynek
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Parameters
t0 = 1e-6       # initial time has to be slightly above 0 to avoid the divergence in the integral

t_final = 0.001    #  arbitrary final time
t = np.linspace(t0, t_final, 1000)# 

# Define the Hubble parameter for the radiation-dominated era
def hubble(time):
    return 1 / (2 * time)

# Numerically integrate H from t0 to t_final
integral_value, error = quad(hubble, t0, t_final)
print("Integral of H from t0 to t_final:", integral_value, "with error estimate:", error)

# Compute n_N at t_final using the integral
n_N_frac =  np.exp(-3 * integral_value)
print("n_N at t =", t_final, "is", n_N_frac)


#analytically

n_N_frac = (t / t0) ** (-3/2)

plt.figure(figsize=(8, 5))
plt.semilogy(t, n_N_frac, label='$n_N(t)$')
plt.xlabel('Time')
plt.ylabel('Number Density $n_N$ (log scale)')
plt.title('Evolution of Number Density/Initial density in Radiation Dominated Era')
plt.legend()
plt.grid(True)
plt.show()