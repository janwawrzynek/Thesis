
from OOP.Particles import Quark, Charged_Lepton, LightNeutrino, SterileNeutrino, Electron, Muon, Tau
from OOP.Universe import Universe
from OOP.Model_Parameters import ModelParameters
from OOP.Rel_degrees_of_freedom import RelDegreesOfFreedom
import numpy as np 
from pathlib import Path
pi = np.pi
exp = np.exp
from scipy.integrate import solve_ivp
from scipy.special import jv, yv
import matplotlib.pyplot as plt
import mpmath
from scipy.integrate import quad
from scipy.special import kn, k1
from scipy.special import zeta

class Simulation:

    def __init__(self, universe: Universe,params:ModelParameters, 
                dof_calculator : RelDegreesOfFreedom):
        
        # ... (rest of your __init__ method)

        #print("Simulation started with", params)
        self.universe = universe
        self.params   = params
        self.dof_calculator = dof_calculator # <-- STORE THE OBJECT
        self.time     = 0
        self.z_list   = np.logspace(-4, 4, num=1000)

        # Pull everything from params instead of hard‑coding:
        self.M_pl = params.M_pl
        self.m_N  = params.m_N
        self.G_F  = params.fermi_constant
        self.f_a  = params.f_a
        self.V_ij = self.params.V_ij
        self.U_matrix = self.params.U_matrix

        

    # NEW MASTER FUNCTION FOR KINEMATIC CHECKS
    def check_kinematics(self, incoming_particle, outgoing_particles):
        """
        Checks if a decay or interaction is kinematically allowed.
        
        Args:
            incoming_particle (Particle): The particle that is decaying.
            outgoing_particles (list[Particle]): A list of the resulting particles.
            
        Raises:
            ValueError: If the sum of outgoing masses is greater than or equal
                        to the incoming particle's mass.
        """

    
        # Get the mass of the incoming particle
        m_incoming = incoming_particle.mass

        # Sum the masses of all outgoing particles
        m_outgoing_total = sum(p.mass for p in outgoing_particles)

        # The core kinematic check
        if m_outgoing_total >= m_incoming:
            raise ValueError(
                "Kinematic condition not satisfied: Sum of outgoing masses "
                f"({m_outgoing_total:.4f} GeV) must be less than the incoming mass "
                f"({m_incoming:.4f} GeV)."
            )
        # If the check passes, the function simply completes.


    # NEW method to get g_star using the provided calculator
    def get_g_star(self, z):
        """
        Calculates the effective number of relativistic degrees of freedom (g_star)
        by calling the provided calculator object.
        """
        return self.dof_calculator.g_star(z)
   
    
    def entropy_density(self, z):
        g_star = self.get_g_star(z)
        return ((2*(pi**2))/(45))*g_star*((m_N/z)**3)
    
    def Hubble(self, z):
        H = (pi/M_pl)*(np.sqrt(self.entropy_density(z)/90))*((self.m_N/z)**2)
        return H
    
    def n_X_eq_bosons(g, mass, z):
        n_eq = zeta(3) * g *( (mass/ z)**3 ) * (1/np.pi **2)

    def n_X_eq_fermions(g, mass, z):
        n_eq = 0.75 * zeta(3) * g *( (mass/ z)**3 ) * (1/np.pi **2)
        return n_eq
    
    def n_X_non_relativistic(g,mass,z):
        return g* ((mass**2)/(2*pi * z))**1.5 * exp(-z)
                              
    def V_ij_for_quark(self, i, j):
        if not isinstance(i, Quark) or not isinstance(j, Quark):
            raise TypeError("Both i and j must be instances of the Quark class.")
        
        # Define quark flavor to index mapping
        up_type_quarks = {'up': 0, 'charm': 1, 'top': 2}  # up, charm, top
        down_type_quarks = {'down': 0, 'strange': 1, 'bottom': 2}  # down, strange, bottom
        
        # Get quark flavors 
        i_flavor = i.flavor.lower()
        j_flavor = j.flavor.lower()
        
        # Determine which quark is up-type and which is down-type
        if i_flavor in up_type_quarks and j_flavor in down_type_quarks:
            i_idx = up_type_quarks[i_flavor]
            j_idx = down_type_quarks[j_flavor]
        elif i_flavor in down_type_quarks and j_flavor in up_type_quarks:
            i_idx = down_type_quarks[i_flavor]
            j_idx = up_type_quarks[j_flavor]
        else:
            raise ValueError("One quark must be up-type (u,c,t) and one must be down-type (d,s,b)")
        
        return abs(self.V_ij[i_idx, j_idx])

        
        
    def lam(self, a,b,c):
        return (a**2) + (b**2) + (c**2) - (2*a*b) - (2*a*c) - (2*b*c)

    
    def Integral(self, x_u,x_d,x_l):
        lower = (x_d + x_l) ** 2
        upper = (1-x_u)**2
        
        def integrand(x):
        # now `x` is a local argument, and x_l, x_d, x_u are closed over
            sqrt_arg = self.lam(x, x_l**2, x_d**2) * self.lam(1, x, x_u**2)

            if sqrt_arg < 0:
                raise ValueError("Square root argument is negative, which is not allowed in this context. It's possible that you have chosen products that aren't possible!")
            
            return (1 / x) * (x - x_l**2 - x_d**2) * (1 + x_u**2 - x) * np.sqrt(sqrt_arg) 
        
        integral, error = quad(integrand, lower, upper)

        return 12 * integral

    def get_U_alpha(self, sterile_neutrino, charged_lepton):
        print("This mixing matrix is a placeholder!!!")
        if not isinstance(sterile_neutrino, SterileNeutrino) and ((isinstance(charged_lepton, Tau)) or isinstance(charged_lepton, Muon) or isinstance(charged_lepton, Electron)):
            raise TypeError("Both sterile_neutrino and charged_lepton must be instances of their respective classes.")
        
        # Matrix columns are electron, muon, tau
        # Matrix rows are n1, n2, n3
        # So if sterile_neutrino is n1 and charged_lepton is electron, then U_alpha[0][0] = 10**(-10)
        # Use self.U_matrix

        # Create mapping:
        flavor_mapping = {
            'n1': 0,
            'n2': 1,
            'n3': 2,
            'electron': 0,
            'muon': 1,
            'tau': 2
        }
        
        # Get indices
        sterile_index = flavor_mapping[sterile_neutrino.flavor]
        charged_lepton_index = flavor_mapping[charged_lepton.__class__.__name__.lower()]
        
        # Return the U_alpha value
        return self.U_matrix[sterile_index][charged_lepton_index]

    
    
    def decay_width_of_sterile_neutrino(self, sterile_neutrino, charged_lepton, u_particle, d_particle):
        # 1) Enforce kinematics before doing anything else:
        self.check_kinematics(
            incoming_particle=sterile_neutrino,
            outgoing_particles=[charged_lepton, u_particle, d_particle]
        )
        # calculate the decay width of a sterile neutrino into a charged lepton and a quark
        U_alpha = self.get_U_alpha(sterile_neutrino, charged_lepton) 
        # generate_U_matrix needs to be defined! Have a library of mixing angles (e.g. e-e mixing angle, e-mu mixing angle, etc.) 

        if isinstance(u_particle, LightNeutrino):
            print("Find better way of implementing this Light Neutrino Mass")
            u_particle.mass = (U_alpha**2) * m_N

        x_l = charged_lepton.mass / self.m_N 
        x_d = d_particle.mass / self.m_N
        x_u = u_particle.mass / self.m_N


        # check type of u_particle (is it a lepton or a quark?)
        if isinstance(u_particle, Charged_Lepton):
            if not isinstance(d_particle, LightNeutrino):
                raise TypeError("For a charged current if u_particle is a Charged_lepton, d_particle must be a LightNeutrino.")
            N_w = 1

        elif isinstance(u_particle, LightNeutrino):
            if not isinstance(d_particle, Charged_Lepton):
                raise TypeError("For Charged current if u_particle is a LightNeutrino, d_particle must be a Charged_Lepton.")
            N_w = 1
        elif isinstance(u_particle, Quark):
            if not isinstance(d_particle, Quark):
                raise TypeError("If u_particle is a Quark, d_particle must also be a Quark.")
            N_c = 3
            N_w = N_c * self.V_ij_for_quark(u_particle, d_particle)**2
        else:
            raise TypeError("u_particle and d_particle must be instances of either the Lepton or Quark class.")

        

        return N_w * (fermi_constant**2) * (self.m_N**5) / (192 * pi**3) * (U_alpha**2) * self.Integral(x_u, x_d, x_l)



if __name__ == "__main__":
    # 1. Define the path to your data file
    # 1. Build a robust, absolute path to the CSV file.
    #    This finds the directory of the current script and looks for the CSV there.
    try:
        script_dir = Path(__file__).parent
        csv_file_path = script_dir / "rel_degree_of_Freedom.csv"
    except NameError:
        # This fallback is useful for interactive environments like Jupyter
        csv_file_path = "rel_degree_of_Freedom.csv"

    
    params     = ModelParameters()    
    # 3. Create the degrees of freedom calculator object
    dof_calculator = RelDegreesOfFreedom(csv_path=csv_file_path, m_N=params.m_N)
  # create the parameters object
    universe = Universe()
    simulation = Simulation(universe, params, dof_calculator)
    
    # Example usage

    e = Electron(0.000511)
    
    #U  = Quark('down') 
    #D  = Quark('up')
    U = LightNeutrino('nu_tau')
    D = Tau(1.777)
    sterile_neutrino = SterileNeutrino('n1')


    # Calculate decay width
    try:
        width = simulation.decay_width_of_sterile_neutrino(sterile_neutrino, e, U, D)
        print(f"Decay width: {width:.6e} GeV")

        # You can also now directly test the g_star function via the simulation object
        g_star_value = simulation.get_g_star(z=1.0)
        print(f"g_star at z=1.0 is: {g_star_value:.2f}")
    except ValueError as ex:
        print("Decay not allowed kinematically:")
        print("  ", ex)
    #decay_width = simulation.decay_width_of_sterile_neutrino(sterile_neutrino, e, U, D)
    #print(f"Decay width of {sterile_neutrino.flavor} into {e.__class__.__name__} and {U.__class__.__name__}  and {D.__class__.__name__}: {decay_width:.6e} GeV")

