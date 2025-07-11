# In OOP/Model.py
import numpy as np
from scipy.integrate import quad
from scipy.special import zeta
from typing import TYPE_CHECKING

# Import your other classes
from .Model_Parameters import ModelParameters
from .Rel_degrees_of_freedom import RelDegreesOfFreedom
# Import particle blueprints
from .Particles import Particle, SterileNeutrino, Charged_Lepton, Electron, Muon, Tau, Quark

# Use TYPE_CHECKING to avoid circular imports at runtime
if TYPE_CHECKING:
    from .Particles import SterileNeutrino, Charged_Lepton

pi = np.pi

class Model:
    """
    Encapsulates a complete physical model, including all parameters,
    constants, and fundamental physical state functions.
    """
    def __init__(self, params: ModelParameters, g_star_csv_path: str):
        print("Initializing physical model...")
        self.params = params
        self.dof_calculator = RelDegreesOfFreedom(
            csv_path=g_star_csv_path, 
            m_N=self.params.m_N
        )

    # --- Particle Factory Methods ---
    def create_sterile_neutrino(self, flavor: str) -> SterileNeutrino:
        """Creates a sterile neutrino with the mass defined in this model."""
        return SterileNeutrino(flavor, mass=self.params.m_N)

    # --- Fundamental State Functions ---
    # These methods describe the properties of the universe under this model.

    def entropy_density(self, z: float) -> float:
        """Calculates the entropy density at a given z."""
        g_star = self.get_g_star(z)
        m_N = self.params.m_N
        return ((2 * (pi**2)) / 45) * g_star * ((m_N / z)**3)
    
    def Hubble(self, z: float) -> float:
        """Calculates the Hubble rate in a radiation-dominated era at a given z."""
        M_pl = self.params.M_pl
        m_N = self.params.m_N
        
        g_star = self.get_g_star(z)
        T = m_N / z
        # Total energy density for radiation
        rho_rad = (pi**2 / 30) * g_star * T**4
        
        H_squared = (8 * pi / (3 * M_pl**2)) * rho_rad
        return np.sqrt(H_squared)

    @staticmethod
    def n_X_eq_bosons(g, mass, z):
        """Calculates equilibrium number density for a bosonic species."""
        return zeta(3) * g * ((mass / z)**3) * (1 / pi**2)

    @staticmethod
    def n_X_eq_fermions(g, mass, z):
        """Calculates equilibrium number density for a fermionic species."""
        return 0.75 * zeta(3) * g * ((mass / z)**3) * (1 / pi**2)

    # --- Helper methods for interactions ---
    
    def get_g_star(self, z: float) -> float:
        return self.dof_calculator.g_star(z)

    def check_kinematics(self, incoming_particle: Particle, outgoing_particles: list[Particle]):
        m_incoming = incoming_particle.mass
        m_outgoing_total = sum(p.mass for p in outgoing_particles)
        if m_outgoing_total >= m_incoming:
            raise ValueError(
                f"Kinematic condition not satisfied: Sum of outgoing masses "
                f"({m_outgoing_total:.4f} GeV) must be less than the incoming mass "
                f"({m_incoming:.4f} GeV)."
            )

    def get_U_alpha(self, sterile_neutrino: 'SterileNeutrino', charged_lepton: 'Charged_Lepton') -> float:
        flavor_mapping = {'n1': 0, 'n2': 1, 'n3': 2, 'electron': 0, 'muon': 1, 'tau': 2}
        sterile_index = flavor_mapping[sterile_neutrino.flavor]
        charged_lepton_index = flavor_mapping[charged_lepton.__class__.__name__.lower()]
        return self.params.U_matrix[sterile_index][charged_lepton_index]

    def V_ij_for_quark(self, i: Quark, j: Quark) -> float:
        up_type = {'up': 0, 'charm': 1, 'top': 2}
        down_type = {'down': 0, 'strange': 1, 'bottom': 2}
        i_flavor, j_flavor = i.flavor.lower(), j.flavor.lower()
        if i_flavor in up_type and j_flavor in down_type:
            i_idx, j_idx = up_type[i_flavor], down_type[j_flavor]
        elif i_flavor in down_type and j_flavor in up_type:
            i_idx, j_idx = down_type[i_flavor], up_type[j_flavor]
        else:
            raise ValueError("One quark must be up-type (u,c,t) and one must be down-type (d,s,b)")
        return abs(self.params.V_ij[i_idx, j_idx])
