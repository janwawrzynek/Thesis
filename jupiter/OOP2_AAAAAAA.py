import numpy as np 
pi = np.pi
exp = np.exp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import jv, yv
import matplotlib.pyplot as plt
import mpmath
from scipy.integrate import quad
from scipy.special import kn, k1
from scipy.special import zeta
# Physical constants (natural units: GeV, s^-1)
M_pl = 2.435e18          # Planck mass [GeV]
fermi_constant = 1.166e-5           # Fermi constant [GeV^-2]
m_N = 0.1 # HNL Mass in GeV
m_a = 1e-6 # ALP mass in GeV
mixing_parameter_squared = 1e-10 # active-sterile neutrino sector mixing
f_a = 1e3 # ALP decay constant scenario 1
mpmath.workdps(50)



class Universe:
    def __init__(self):
        self.Particles = []
    
    def add_Particle(self, Particle):
        self.Particles.append(Particle)
    
    def remove_Particle(self, Particle):
        if Particle in self.Particles: # n time :(
            self.Particles.remove(Particle)

    def __str__(self):
        return f"Universe: {self.alpha}, Particles: {len(self.Particles)}"
    

class Particle:
    def __init__(self, mass):
        self.mass = mass
    
    def __str__(self):
        return f"Particle: {self.alpha, self.mass}" 
    
'''
class Quark(Particle):
    def __init__(self, flavor):
        super().__init__()
        
        if flavor not in ['up', 'down', 'charm', 'strange', 'top', 'bottom']:
            raise ValueError("Invalid quark flavor. Must be one of: up, down, charm, strange, top, bottom.")
        self.flavor = flavor

        if flavor == 'up':
            self.mass = 0.005
        elif flavor == 'down':
            self.mass = 0.003
        elif flavor == 'charm':
            self.mass = 1.3
        elif flavor == 'strange':
            self.mass = 0.1
        elif flavor == 'top':
            self.mass = 174 
        elif flavor == 'bottom':
            self.mass = 4.5


    def __str__(self):
        return f"Quark: {self.flavor}, Mass: {self.mass}"
    
'''

class Quark(Particle):
    def __init__(self, flavor):
        # 1) Validate flavor
        if flavor not in ['up','down','charm','strange','top','bottom']:
            raise ValueError("…")

        # 2) Look up the correct numerical mass
        mass_map = {
            'up':      0.005,
            'down':    0.003,
            'charm':   1.3,
            'strange': 0.1,
            'top':     174,
            'bottom':  4.5
        }
        mass_value = mass_map[flavor]

        # 3) Now initialize the parent with that mass
        super().__init__(mass_value)

        # 4) Store your flavor
        self.flavor = flavor

class Lepton(Particle):
    def __init__(self, mass):
        Particle.__init__(self, mass) 

class Charged_Lepton(Lepton):
    def __init__(self, mass):
        super().__init__(mass)

        

class Neutrino(Lepton):
    def __init__(self):
        super().__init__(mass=0.0)

class LightNeutrino(Neutrino):
    def __init__(self, flavor):
        super().__init__()
        
        if flavor not in ['nu_e', 'nu_mu', 'nu_tau']:
            raise ValueError("Invalid neutrino flavor. Must be one of: nu_e, nu_mu, nu_tau.")
        
        self.flavor = flavor
        print("CHECK THE LIGHT NEUTRINO COMPOSITION. SETTING THE MASS TO 0.001 FOR NOW")
        self.mass = 0.001  # Light neutrino mass in GeV

    def __str__(self):
        return f"Neutrino: {self.flavor}, Mass: {self.mass}"

class SterileNeutrino(Neutrino):
    def __init__(self, flavor):
        super().__init__()
        
        if flavor not in ['n1', 'n2', 'n3']:
            raise ValueError("Invalid sterile neutrino flavor. Must be one of: n1, n2, n3.")

        self.flavor = flavor
        print("CHECK THE STERILE NEUTRINO COMPOSITION. SETTING THE MASS TO 0.1 FOR NOW")
        self.mass = 0.1


    def __str__(self):
        return f"Quark: {self.flavor}, Mass: {self.mass}"

class Muon(Charged_Lepton):
    def __init__(self, mass):
        super().__init__(mass)
        print("CHECK THE MUON COMPOSITION. SETTING THE MASS TO 0.1056583745 FOR NOW")
        self.mass = 0.106  # Muon mass in GeV

class Electron(Charged_Lepton):
    def __init__(self, mass):
        super().__init__(mass)
        print("CHECK THE ELECTRON COMPOSITION. SETTING THE MASS TO 0.000511 FOR NOW")
        self.mass = 0.0005  # Electron mass in GeV

class Tau(Charged_Lepton):
    def __init__(self, mass):
        super().__init__(mass)
        print("CHECK THE TAU COMPOSITION. SETTING THE MASS TO 1.77686 FOR NOW")
        self.mass = 1.78  # Tau mass in GeV



    

class Simulation:
    def __init__(self, universe):
        print("Simulation started")
        self.universe = universe
        self.time = 0
        self.z_list = np.logspace(-4, 4, num=1000)
        self.M_pl = 2.435e18  # Planck mass in GeV
        self.m_N = 0.1 # HNL mass in GeV
        self.V_ij = np.array([[0.974, 0.225, 0.004], [0.225, 0.973, 0.041], [0.009, 0.04, 0.999]])
        self.U_matrix = np.array([[10**(-10), 10**(-10), 0], [0, 10**(-10), 0], [0, 0, 10**(-10)]])

    
    def get_g_star(self, z):
        if z > 0.1/173.3 and z < 0.1/125.6 :
            g_star = 106.75 

        elif z > 0.1/125.6 and z < 0.1/91.2 :

            g_star = 96.25
        elif z> 0.1/91.2 and z <  0.1/80.4:
            g_star = 95.25
        elif z> 0.1/80.4 and z < 0.1/4.19:
            g_star = 92.25
        elif z> 0.1/4.19 and z < 0.1/1.777 :
            g_star = 86.25
        elif z>0.1/1.777 and  z < 0.1/1.29:
            g_star = 75.75
        elif z> 0.1/1.29 and z < 0.1/0.214 :
            g_star = 72.25
        elif z> 0.1/0.214 and z < 0.1/0.1396:
            g_star = 61.75
        elif z> 0.1/0.1396 and z < 0.1/0.135 :
            g_star = 17.25
        elif z> 0.1/0.135 and z < 0.1/0.1057 :
            g_star = 15.25
        elif z> 0.1/0.1057 and z < 0.1/0.01:
            g_star = 14.25
        elif z>0.1/800e-6 and z < 0.1/511e-3:
            g_star = 10.75
        elif z < 800e-6:
            g_star = 3.909
        
        return g_star
    
    
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
        # calculate the decay width of a sterile neutrino into a charged lepton and a quark
        U_alpha = self.get_U_alpha(sterile_neutrino, charged_lepton) 
        # generate_U_matrix needs to be defined! Have a library of mixing angles (e.g. e-e mixing angle, e-mu mixing angle, etc.) 

        if isinstance(u_particle, LightNeutrino):
            print("Find better way of implementing this Light Neutrino Mass")
            u_particle.mass = (U_alpha**2) * m_N

        x_l = charged_lepton.mass / self.m_N 
        x_d = d_particle.mass / self.m_N
        x_u = u_particle.mass / self.m_N

        # Kinematic check:
        if x_u + x_d + x_l >= 1:
            raise ValueError("Kinematic condition not satisfied: x_u + x_d + x_l must be less than 1.")



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
    universe = Universe()
    simulation = Simulation(universe)

    # Example usage

    e = Electron(0.000511)
    
    #U  = Quark('down') 
    #D  = Quark('up')
    U = LightNeutrino('nu_e')
    D = Electron(0.000511)  
    sterile_neutrino = SterileNeutrino('n1')

    # Calculate decay width
    decay_width = simulation.decay_width_of_sterile_neutrino(sterile_neutrino, e, U, D)
    print(f"Decay width of {sterile_neutrino.flavor} into {e.__class__.__name__} and {U.__class__.__name__}  and {D.__class__.__name__}: {decay_width:.6e} GeV")



