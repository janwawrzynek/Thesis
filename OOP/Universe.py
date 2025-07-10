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
    
