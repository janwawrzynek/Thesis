
import numpy as np
class ModelParameters:
    def __init__(self,
                 M_pl=2.435e18,
                 fermi_constant=1.166e-5,
                 m_N=0.1,
                 mixing_parameter_squared=1e-10,
                 f_a=1e3,

                 # ← new defaults for your matrices:
                 V_ij=None,
                 U_matrix=None):
        self.M_pl = M_pl
        self.fermi_constant = fermi_constant
        self.m_N = m_N
        self.mixing_parameter_squared = mixing_parameter_squared
        self.f_a = f_a

        # CKM matrix for quarks:
        self.V_ij = (np.array(V_ij)
                     if V_ij is not None
                     else np.array([
                         [0.974, 0.225, 0.004],
                         [0.225, 0.973, 0.041],
                         [0.009, 0.040, 0.999],
                     ]))

        # active–sterile mixing:
        self.U_matrix = (np.array(U_matrix)
                         if U_matrix is not None
                         else np.array([
                             [1e-10, 1e-10, 0],
                             [0,     1e-10, 0],
                             [0,     0,     1e-10],
                         ]))

    def __repr__(self):
        return (f"ModelParameters(m_N={self.m_N}, G_F={self.fermi_constant}, "
                f"V_ij=\n{self.V_ij},\nU_matrix=\n{self.U_matrix})")
