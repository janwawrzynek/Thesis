# In OOP/Particles.py
import numpy as np
from typing import List, TYPE_CHECKING
from scipy.integrate import quad

# This now imports the Model class for type hinting
if TYPE_CHECKING:
    from OOP.Model import Model

pi = np.pi

# --- Base Particle Class ---
class Particle:
    def __init__(self, mass):
        self.mass = mass
    def __str__(self):
        return f"{self.__class__.__name__} (Mass: {self.mass} GeV)"

# --- Lepton Base Classes ---
class Lepton(Particle): pass
class Charged_Lepton(Lepton): pass
class Neutrino(Lepton): pass

# --- Specific Charged Leptons ---
class Electron(Charged_Lepton):
    def __init__(self):
        super().__init__(mass=0.000511)

class Muon(Charged_Lepton):
    def __init__(self):
        super().__init__(mass=0.105658)

class Tau(Charged_Lepton):
    def __init__(self):
        super().__init__(mass=1.77686)

# --- Specific Neutrinos ---
class LightNeutrino(Neutrino):
    def __init__(self, flavor):
        super().__init__(mass=0.001) 
        if flavor not in ['nu_e', 'nu_mu', 'nu_tau']:
            raise ValueError("Invalid light neutrino flavor.")
        self.flavor = flavor

class SterileNeutrino(Neutrino):
    def __init__(self, flavor, mass=0.1):
        super().__init__(mass)
        if flavor not in ['n1', 'n2', 'n3']:
            raise ValueError("Invalid sterile neutrino flavor.")
        self.flavor = flavor

    def __str__(self):
        return f"SterileNeutrino: {self.flavor}, Mass: {self.mass}"

    # --- Decay Calculation Methods ---
    # These methods are now part of the particle's intrinsic definition.

    @staticmethod
    def lam(a, b, c):
        """Källén function for kinematic calculations. Static as it's a pure function."""
        return (a**2) + (b**2) + (c**2) - (2*a*b) - (2*a*c) - (2*b*c)

    def _integral(self, x_u, x_d, x_l):
        """Internal phase space integral for 3-body decay."""
        lower = (x_d + x_l) ** 2
        upper = (1 - x_u) ** 2
        
        # Define the integrand within the method's scope
        def integrand(x):
            # It can now call self.lam directly
            sqrt_arg = self.lam(x, x_l**2, x_d**2) * self.lam(1, x, x_u**2)
            if sqrt_arg < 0: return 0.0
            return (1 / x) * (x - x_l**2 - x_d**2) * (1 + x_u**2 - x) * np.sqrt(sqrt_arg)
            
        integral, _ = quad(integrand, lower, upper)
        return 12 * integral

    def calculate_decay_width(self, outgoing_particles: List[Particle], model: 'Model') -> float:
        """
        Calculates the decay width using the provided physical model for context.
        """
        # --- 1. Use the model's rules to check kinematics FIRST ---
        model.check_kinematics(self, outgoing_particles)

        # --- 2. Unpack particles and validate channel ---
        if len(outgoing_particles) != 3:
            raise ValueError("This decay channel is defined for 3 outgoing particles.")
        charged_lepton, u_particle, d_particle = outgoing_particles

        # --- 3. Get parameters from the model ---
        U_alpha = model.get_U_alpha(self, charged_lepton)
        m_N = model.params.m_N
        G_F = model.params.fermi_constant
        
        if isinstance(u_particle, LightNeutrino):
            u_particle.mass = (U_alpha**2) * m_N

        x_l = charged_lepton.mass / m_N
        x_d = d_particle.mass / m_N
        x_u = u_particle.mass / m_N

        # --- 4. Determine the interaction factor N_w using model helpers ---
        if isinstance(u_particle, Charged_Lepton):
            if not isinstance(d_particle, LightNeutrino):
                raise TypeError("For charged current, if u_particle is a ChargedLepton, d_particle must be a LightNeutrino.")
            N_w = 1
        elif isinstance(u_particle, LightNeutrino):
            if not isinstance(d_particle, Charged_Lepton):
                raise TypeError("For charged current, if u_particle is a LightNeutrino, d_particle must be a ChargedLepton.")
            N_w = 1
        elif isinstance(u_particle, Quark):
            if not isinstance(d_particle, Quark):
                raise TypeError("If u_particle is a Quark, d_particle must also be a Quark.")
            N_c = 3
            N_w = N_c * model.V_ij_for_quark(u_particle, d_particle)**2
        else:
            raise TypeError("u_particle and d_particle must be instances of Lepton or Quark.")

        # --- 5. Calculate the final width using its OWN integral method ---
        integral_val = self._integral(x_u, x_d, x_l)
        width = N_w * (G_F**2) * (m_N**5) / (192 * pi**3) * (U_alpha**2) * integral_val
        return width

# --- Quark Class ---
class Quark(Particle):
    def __init__(self, flavor):
        if flavor not in ['up','down','charm','strange','top','bottom']:
            raise ValueError("Invalid quark flavor.")
        mass_map = {'up':0.005, 'down':0.003, 'charm':1.3, 'strange':0.1, 'top':174, 'bottom':4.5}
        super().__init__(mass_map[flavor])
        self.flavor = flavor
