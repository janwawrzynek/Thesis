# In OOP/Model_Parameters.py
import numpy as np

class ModelParameters:
    """
    A container for all the fundamental constants and user-definable
    parameters of a given physical model.
    """
    def __init__(self):
        # --- Physical Constants ---
        self.M_pl = 2.435e18          # Planck mass [GeV]
        self.fermi_constant = 1.166e-5   # Fermi constant [GeV^-2]
        
        # ADDED: The sine-squared of the Weinberg angle, needed for NC calculations.
        self.sin_theta_w_sq = 0.2312 

        # --- Model-specific Parameters (can be changed after initialization) ---
        self.m_N = 0.1 # HNL Mass in GeV (default)
        self.f_a = 1e3 # ALP decay constant (example)

        # --- Mixing Matrices (Placeholders) ---
        # CKM Matrix
        self.V_ij = np.array([[0.974, 0.225, 0.004],
                              [0.225, 0.973, 0.041],
                              [0.009, 0.040, 0.999]])
        
        # PMNS-like mixing matrix for HNLs
        self.U_matrix = np.array([[1e-5, 1e-5, 0],
                                  [0, 1e-5, 0],
                                  [0, 0, 1e-5]])
