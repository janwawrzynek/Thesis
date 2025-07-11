import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import jv, yv
import matplotlib.pyplot as plt
import mpmath
from scipy.integrate import quad
from scipy.special import kn, k1
# Physical constants (natural units: GeV, s^-1)
M_pl = 2.435e18          # Planck mass [GeV]
fermi_constant = 1.166e-5           # Fermi constant [GeV^-2]
m_N = 0.1 # HNL Mass in GeV
m_a = 1e-6 # ALP mass in GeV
mixing_parameter_squared = 1e-10 # active-sterile neutrino sector mixing
f_a = 1e6 # ALP decay constant scenario 1
#f_M = 0.1 # in GeV
m_nu = mixing_parameter_squared*m_N
z_init = 1e-4
z_final = 1e2
mass_temp_ratio=  np.linspace(z_init,z_final, 100)
def decay_width_N_av(f_a,m_N,mixing_parameter_squared, m_a):
    decay_width_N_av = mixing_parameter_squared*(m_N**3)*(1/(4*np.pi*(f_a**2)))*(np.sqrt(1+ (m_a/m_N)**2))*((1-(m_a/m_N)**2)**1.5)
    return decay_width_N_av
def decay_width_a_vv(f_a, m_N, mixing_parameter_squared, m_a, m_v):
    decay_width_a_vv = (1/(2* np.pi * f_a**2))*m_N**2 * m_a *( mixing_parameter_squared**2) * (np.sqrt(1- (4*m_nu**2)/m_a**2)) * (1- (2*m_nu**2)/m_a**2)
    #decay_width_a_vv = ((m_N**2 )* m_a * (mixing_parameter_squared**2))/(2* np.pi * f_a**2)
    return decay_width_a_vv
gamma_N_av_n_N_eq_ratio_list = [] 
gamma_a_vv_n_N_eq_ratio_list = []
for z in mass_temp_ratio:
    gamma_N_av_n_N_eq_ratio_value =  decay_width_N_av(1e3,m_N,mixing_parameter_squared, m_a)*(kn(1,z)/kn(2,z))
    gamma_N_av_n_N_eq_ratio_list.append(gamma_N_av_n_N_eq_ratio_value)
    gamma_a_vv_n_N_eq_ratio_value = decay_width_a_vv(f_a, m_N, mixing_parameter_squared, m_a, m_nu)*(kn(1,z)/kn(2,z))
    gamma_a_vv_n_N_eq_ratio_list.append(gamma_a_vv_n_N_eq_ratio_value)
plt.loglog(mass_temp_ratio, gamma_N_av_n_N_eq_ratio_list, label = r'$\gamma_{N \rightarrow a \nu }/ n_{N}^{eq}$')
plt.loglog(mass_temp_ratio, gamma_a_vv_n_N_eq_ratio_list,label = r'$\gamma_{a \rightarrow \nu \nu }/ n_{a}^{eq}$' )
plt.xlim(1e-4, 1e2)
plt.title(r'$\gamma_{X \rightarrow Y }/ n_{X}^{eq}$')
plt.xlabel("z = m_N /T")
plt.ylabel(r'$\gamma_{X \rightarrow Y }/ n_{X}^{eq}$')
plt.legend()
plt.show()
