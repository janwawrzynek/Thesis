# In OOP/Simulation.py
import numpy as np
from pathlib import Path

# Import the new Model class and the particle blueprints
from .Model import Model
from .Model_Parameters import ModelParameters
from .Particles import Electron, Tau, LightNeutrino, Muon

pi = np.pi
exp = np.exp

class Simulation:
    """
    The simulation engine. It takes a physical model and calculates
    how a system governed by that model evolves over time.
    """
    def __init__(self, model: Model):
        """
        Initializes the simulation engine with a specific physical model.
        
        Args:
            model (Model): The self-contained physical model to simulate.
        """
        print("Initializing simulation engine...")
        self.model = model
        self.z_list = np.logspace(-4, 4, num=1000)

    def run(self):
        """
        This method will eventually contain the logic to solve the
        Boltzmann equations over the z_list, using the model's
        properties like Hubble rate and decay widths.
        """
        print("\nSimulation run method called. (Boltzmann solver would go here)")
        # The following diagnostic loop has been commented out to prevent printing the z, H, and s values.
        # for z in self.z_list[:5]: # Just showing first 5 steps
        #     H = self.model.Hubble(z)
        #     s = self.model.entropy_density(z)
        #     print(f"z = {z:.2f}, H = {H:.4e} GeV, s = {s:.4e} GeV^3")
        

# --- Main Execution Block ---
if __name__ == "__main__":
    
    # --- 1. Define the Model ---
    # First, set up the parameters for the specific model you want to test.
    params = ModelParameters()
    params.m_N = 0.2 # Example: testing a 0.5 GeV sterile neutrino
    
    # Define the path to the g* data file
    try:
        script_dir = Path(__file__).parent
        # Go up one level from OOP to the Thesis directory to find the CSV
        csv_path = script_dir.parent / "rel_degree_of_Freedom.csv"
        if not csv_path.is_file():
             # Fallback for different structures
             csv_path = script_dir / "rel_degree_of_Freedom.csv"
    except NameError:
        csv_path = "rel_degree_of_Freedom.csv"
        
    # Now, create the complete Model object
    my_model = Model(params=params, g_star_csv_path=str(csv_path))

    # --- 2. Set up the Simulation ---
    # The simulation engine just takes the model as a whole.
    simulation = Simulation(model=my_model)
    # Run the simulation
    simulation.run()

    # --- 3. Run a Specific Calculation (like a decay width) ---
    # Use the model's "factory" to create a sterile neutrino instance
    # that is consistent with the model's parameters (e.g., mass).
    sterile_neutrino = my_model.create_sterile_neutrino('n1')
    
    # Create the other particles involved in the decay
    e = Electron()
    U = LightNeutrino('nu_mu')
    D = Muon() # Example decay to a muon
    
    print(f"\nCalculating decay for a {sterile_neutrino.mass} GeV Sterile Neutrino.")

    try:
        # Ask the particle to calculate its decay width, providing the model as context.
        outgoing_particles = [e, U, D]
        width = sterile_neutrino.calculate_decay_width(outgoing_particles, my_model)

        print(f"Decay width of {sterile_neutrino.flavor}: {width:.6e} GeV")

    except ValueError as ex:
        print("Decay not allowed:")
        print(f"   {ex}")
