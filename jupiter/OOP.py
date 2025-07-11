class Universe:
    def __init__(self):
        self.temperature = 
        self.z = 
        self.particles = []
    
    def add_particle(self, particle):
        self.particles.append(particle)
    
    def remove_particle(self, particle):
        if particle in self.particles: # n time :(
            self.particles.remove(particle)

    def __str__(self):
        return f"Universe: {self.name}, Particles: {len(self.particles)}"



class Particle:
    def __init__(self, name, mass, decay_const, charge=None):
        self.name = name
        self.mass = mass
        self.decay_const = decay_const
        self.charge = charge
    
    def __str__(self):
        return f"ParticleA: {self.name}, Mass: {self.mass}, Charge: {self.charge}"

    # def Y(self):
    #     y_N = self.universe.get_Y_value("N", old_Y) # ask the universe for Y value of particle type N, given old_Y


class Lepton(Particle):
    def __init__(self, name, mass, decay_const, charge):
        super().__init__(name, mass, decay_const)
        self.charge = charge
    
    def __str__(self):
        return f"Lepton: {self.name}, Mass: {self.mass}, Charge: {self.charge}"
    
    # gamma func for lepton with ll
    # gamma func for lepton with pp



    def gamma(self, other_particle):
        if isinstance(other_particle, Lepton):
            # Calculate gamma for lepton with lepton
            return self.decay_const * other_particle.mass / self.mass # l_mass * (l_mass / l_mass)
        elif isinstance(other_particle, Proton):
            # Calculate gamma for lepton with other particle
            # l_mass * (p_mass / l_mass) + (l_charge * p_charge)
            return self.decay_const * other_particle.mass / self.mass + (self.charge * other_particle.charge / self.charge)
        else:
            raise ValueError("Invalid particle type for gamma calculation")

class Proton(Particle):
    def __init__(self, name, mass, decay_const, charge):
        super().__init__(name, mass, decay_const)
        self.charge = charge
    
    def __str__(self):
        return f"Proton: {self.name}, Mass: {self.mass}, Charge: {self.charge}"
class Meson(Particle):
    def __init__(self, name, mass, decay_const, charge, coupling):
        super().__init__(name, mass, decay_const, charge)
        self.coupling = coupling
    
    def __str__(self):
        return f"Meson: {self.name}, Mass: {self.mass}, Charge: {self.charge}, Coupling: {self.coupling}"
     

# class Electron(Lepton):
#     def __init__(self, name, mass, decay_const, charge):
#         super().__init__(name, mass, decay_const, charge)
    
#     def __str__(self):
#         return f"Electron: {self.name}, Mass: {self.mass}, Charge: {self.charge}"
# class Muon(Lepton):
#     def __init__(self, name, mass, decay_const, charge):
#         super().__init__(name, mass, decay_const, charge)
    
#     def __str__(self):
#         return f"Muon: {self.name}, Mass: {self.mass}, Charge: {self.charge}"
# class Tau(Lepton):
#     def __init__(self, name, mass, decay_const, charge):
#         super().__init__(name, mass, decay_const, charge)
    
#     def __str__(self):
#         return f"Tau: {self.name}, Mass: {self.mass}, Charge: {self.charge}"

# class Neutrino(Lepton):
#     def __init__(self, name, mass, decay_const, charge):
#         super().__init__(name, mass, decay_const, charge)
#         self.t = Tau(name, decay_const, charge)
#         self.m = Muon(name, decay_const, charge)
#         self.e = Electron(name, decay_const, charge)
    
#     def interact(self, other_particle):
#         random_value = random.random()
#         if random_value < 0.3:
#             return self.t.gamma(other_particle)


# Worth defining particles as separate classes
# This allows for bundling of attributes and methods specific to each particle type
# And it allows for interactions of particle types with particles of other types
# And that allows for bundling of very specific functions to a particle
# e.g. gamma_N_to_SM
class Neutron(Particle):
    def __init__(self, name, mass, decay_const, charge):
        super().__init__(name, mass, decay_const, charge)
    def __str__(self):
        return f"Neutron: {self.name}, Mass: {self.mass}, Charge"
        return self.charge
    def gamma(self, other_particle):
        if isinstance(other_particle, Proton):
            # Calculate gamma for neutron with proton
            return self.decay_const * other_particle.mass / self.mass + (self.charge * other_particle.charge / self.charge)
        elif isinstance(other_particle, Neutron):
            # Calculate gamma for neutron with neutron
            return self.decay_const * other_particle.mass / self.mass
        else:
            raise ValueError("Invalid particle type for gamma calculation")