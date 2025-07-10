

class Particle:
    def __init__(self, mass):
        self.mass = mass
    
    def __str__(self):
        return f"Particle: {self.alpha, self.mass}" 
    

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
        #print("CHECK THE LIGHT NEUTRINO COMPOSITION. SETTING THE MASS TO 0.001 FOR NOW")
        self.mass = 0.001  # Light neutrino mass in GeV

    def __str__(self):
        return f"Neutrino: {self.flavor}, Mass: {self.mass}"

class SterileNeutrino(Neutrino):
    def __init__(self, flavor):
        super().__init__()
        
        if flavor not in ['n1', 'n2', 'n3']:
            raise ValueError("Invalid sterile neutrino flavor. Must be one of: n1, n2, n3.")

        self.flavor = flavor
        #print("CHECK THE STERILE NEUTRINO COMPOSITION. SETTING THE MASS TO 0.1 FOR NOW")
        self.mass = 0.1


    def __str__(self):
        return f"Quark: {self.flavor}, Mass: {self.mass}"

class Muon(Charged_Lepton):
    def __init__(self, mass):
        super().__init__(mass)
        #print("CHECK THE MUON COMPOSITION. SETTING THE MASS TO 0.1056583745 FOR NOW")
        self.mass = 0.106  # Muon mass in GeV

class Electron(Charged_Lepton):
    def __init__(self, mass):
        super().__init__(mass)
        #print("CHECK THE ELECTRON COMPOSITION. SETTING THE MASS TO 0.000511 FOR NOW")
        self.mass = 0.0005  # Electron mass in GeV

class Tau(Charged_Lepton):
    def __init__(self, mass):
        super().__init__(mass)
        #print("CHECK THE TAU COMPOSITION. SETTING THE MASS TO 1.77686 FOR NOW")
        self.mass = 1.78  # Tau mass in GeV

