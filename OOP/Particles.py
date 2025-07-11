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

# --- Quark Class ---
class Quark(Particle):
    def __init__(self, flavor):
        if flavor not in ['up','down','charm','strange','top','bottom']:
            raise ValueError("Invalid quark flavor.")
        mass_map = {'up':0.005, 'down':0.003, 'charm':1.3, 'strange':0.1, 'top':174, 'bottom':4.5}
        super().__init__(mass_map[flavor])
        self.flavor = flavor
        self.is_up_type = flavor in ['up', 'charm', 'top']

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
    @staticmethod
    def lam(a, b, c):
        """Källén function for kinematic calculations."""
        return (a**2) + (b**2) + (c**2) - (2*a*b) - (2*a*c) - (2*b*c)

    def _integral_cc(self, x_u, x_d, x_l):
        """Internal phase space integral for 3-body CHARGED current decay."""
        lower = (x_d + x_l) ** 2
        upper = (1 - x_u) ** 2
        def integrand(x):
            sqrt_arg = self.lam(x, x_l**2, x_d**2) * self.lam(1, x, x_u**2)
            if sqrt_arg < 0: return 0.0
            return (1 / x) * (x - x_l**2 - x_d**2) * (1 + x_u**2 - x) * np.sqrt(sqrt_arg)
        integral, _ = quad(integrand, lower, upper)
        return 12 * integral

    def _get_nc_coeffs(self, fermion: Particle, is_interference: bool, model: 'Model'):
        """
        Calculates the C1 and C2 coefficients from Table 3 of the paper.
        """
        s2w = model.params.sin_theta_w_sq # sin^2(theta_W)
        if isinstance(fermion, Quark):
            if fermion.is_up_type:
                C1 = 0.25 * (1 - (8/3)*s2w + (32/9)*s2w**2)
                C2 = (1/3)*s2w * ((4/3)*s2w - 1)
            else: # down-type quark
                C1 = 0.25 * (1 - (4/3)*s2w + (8/9)*s2w**2)
                C2 = (1/6)*s2w * ((4/3)*s2w - 1)
        elif isinstance(fermion, Charged_Lepton):
            if is_interference:
                C1 = 0.25 * (1 + 4*s2w + 8*s2w**2)
                C2 = 0.5 * s2w * (2*s2w + 1)
            else:
                C1 = 0.25 * (1 - 4*s2w + 8*s2w**2)
                C2 = 0.5 * s2w * (2*s2w - 1)
        else:
            raise TypeError(f"Cannot calculate NC coeffs for particle type {fermion.__class__.__name__}")
        return C1, C2
    
    @staticmethod
    def _L_function(x):
        """The L(x) function from Eq. 3.4 of the paper, made numerically stable."""
        if x >= 0.5: # Kinematically forbidden
            return 0
        if x == 0:
            return 0
        
        sqrt_val = np.sqrt(1 - 4*x**2)
        numerator = 1 - 3*x**2 - (1 - x**2) * sqrt_val
        denominator = x**2 * (1 + sqrt_val)

        # Prevent log of zero or negative from floating point errors
        if numerator <= 0 or denominator <= 0:
            return 0
            
        arg = numerator / denominator
        if arg <= 0:
            return 0
            
        return np.log(arg)

    def calculate_decay_width(self, outgoing_particles: List[Particle], model: 'Model') -> float:
        """
        Calculates the decay width for various 3-body channels using the provided model.
        This function now distinguishes between charged and neutral current processes.
        """
        model.check_kinematics(self, outgoing_particles)
        if len(outgoing_particles) != 3:
            raise ValueError("This decay channel is defined for 3 outgoing particles.")
        
        p1, p2, p3 = outgoing_particles
        m_N = model.params.m_N
        G_F = model.params.fermi_constant
        
        # Determine the primary lepton flavor for U_alpha mixing.
        # This is the charged lepton that couples to the incoming HNL at the vertex.
        primary_lepton = None
        if isinstance(p1, Charged_Lepton): # CC decay: N -> l_alpha + ...
            primary_lepton = p1
        elif isinstance(p1, LightNeutrino): # NC decay: N -> nu_alpha + ...
            # The mixing is with the outgoing active neutrino's flavor.
            flavor_map = {'nu_e': Electron(), 'nu_mu': Muon(), 'nu_tau': Tau()}
            primary_lepton = flavor_map.get(p1.flavor)
        
        if primary_lepton is None:
             raise TypeError("Could not determine primary lepton flavor for U_alpha mixing.")

        U_alpha = model.get_U_alpha(self, primary_lepton)

        # --- Case 1: Charged Current Decay (N -> l_alpha + X + Y) ---
        if isinstance(p1, Charged_Lepton):
            charged_lepton, u_particle, d_particle = p1, p2, p3
            x_l = charged_lepton.mass / m_N
            x_d = d_particle.mass / m_N
            x_u = u_particle.mass / m_N
            
            if isinstance(u_particle, LightNeutrino) and isinstance(d_particle, Charged_Lepton):
                lepton_doublets = {'nu_e': Electron, 'nu_mu': Muon, 'nu_tau': Tau}
                expected_partner = lepton_doublets.get(u_particle.flavor)
                if not isinstance(d_particle, expected_partner):
                    raise TypeError(f"Lepton family violation: {u_particle.flavor} must pair with {expected_partner.__name__}.")
                if charged_lepton.__class__ == d_particle.__class__:
                    raise ValueError("This is an interference case (alpha=beta) and must be calculated as a neutral current decay.")
                N_w = 1
            elif isinstance(u_particle, Quark) and isinstance(d_particle, Quark):
                N_c = 3
                N_w = N_c * model.V_ij_for_quark(u_particle, d_particle)**2
            else:
                raise TypeError("Invalid particle combination for charged current decay.")
            
            integral_val = self._integral_cc(x_u, x_d, x_l)
            return N_w * (G_F**2) * (m_N**5) / (192 * pi**3) * (U_alpha**2) * integral_val

        # --- Case 2: Neutral Current Decay (N -> nu_alpha + f + f_bar) ---
        elif isinstance(p1, LightNeutrino):
            neutrino, fermion1, fermion2 = p1, p2, p3
            
            if isinstance(fermion1, LightNeutrino) and isinstance(fermion2, LightNeutrino):
                delta_ab = 1.0 if neutrino.flavor == fermion1.flavor else 0.0
                return (1 + delta_ab) * (G_F**2 * m_N**5) / (768 * pi**3) * abs(U_alpha)**2
            
            if fermion1.__class__ != fermion2.__class__:
                raise TypeError("Neutral current decay products must be a particle-antiparticle pair.")
            
            x = fermion1.mass / m_N
            if x >= 0.5: return 0.0 # Kinematically forbidden
            
            sqrt_term = np.sqrt(1 - 4*x**2)
            L_val = self._L_function(x)
            
            lepton_flavor_map = {Electron: 'nu_e', Muon: 'nu_mu', Tau: 'nu_tau'}
            is_interference = lepton_flavor_map.get(fermion1.__class__) == neutrino.flavor

            C1, C2 = self._get_nc_coeffs(fermion1, is_interference, model)
            
            term1 = C1 * ((1 - 14*x**2 - 2*x**4 - 12*x**6) * sqrt_term + 12*x**4 * (x**4 - 1) * L_val)
            term2 = 4 * C2 * (x**2 * (2 + 10*x**2 - 12*x**4) * sqrt_term + 6*x**4 * (1 - 2*x**2 + 2*x**4) * L_val)
            
            N_z = 3 if isinstance(fermion1, Quark) else 1
            
            return N_z * (G_F**2 * m_N**5) / (192 * pi**3) * abs(U_alpha)**2 * (term1 + term2)

        else:
            raise TypeError("Decay must be initiated by either a charged lepton or a neutrino.")
