# In OOP/Simulation.py
import numpy as np
from pathlib import Path

# Import the new Model class and the particle blueprints
from .Model import Model
from .Model_Parameters import ModelParameters
from .Particles import Electron, Tau, LightNeutrino, Muon, Quark

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

# --- Main Execution Block ---
if __name__ == "__main__":
    
    # --- 1. Define the Model ---
    params = ModelParameters()
    params.m_N = 0.5 # Example: testing a 0.5 GeV sterile neutrino
    
    try:
        script_dir = Path(__file__).parent
        csv_path = script_dir.parent / "rel_degree_of_Freedom.csv"
        if not csv_path.is_file():
             csv_path = script_dir / "rel_degree_of_Freedom.csv"
    except NameError:
        csv_path = "rel_degree_of_Freedom.csv"
        
    my_model = Model(params=params, g_star_csv_path=str(csv_path))

    # --- 2. Set up the Simulation ---
    simulation = Simulation(model=my_model)
    simulation.run()

    # --- 3. Run Specific Calculations ---
    sterile_neutrino = my_model.create_sterile_neutrino('n1')
    print(f"\n--- Calculating decays for a {sterile_neutrino.mass} GeV Sterile Neutrino ---")

    # --- CALCULATION A: Charged Current Decay (alpha != beta) ---
    print("\n[A] Charged Current Decay: N -> e- + nu_mu + mu+")
    try:
        # For N -> l_alpha + X + Y, the order is [l_alpha, X, Y]
        outgoing_cc = [Electron(), LightNeutrino('nu_mu'), Muon()]
        width_cc = sterile_neutrino.calculate_decay_width(outgoing_cc, my_model)
        print(f"  > Decay width: {width_cc:.6e} GeV")
    except (ValueError, TypeError) as ex:
        print(f"  > Calculation failed: {ex}")

    # --- CALCULATION B: Neutral Current + Interference Decay (alpha = beta) ---
    print("\n[B] Neutral Current (Interference) Decay: N -> nu_e + e+ + e-")
    try:
        # For N -> nu + f + f_bar, the order is [nu, f, f_bar]
        outgoing_nc = [LightNeutrino('nu_e'), Electron(), Electron()]
        width_nc = sterile_neutrino.calculate_decay_width(outgoing_nc, my_model)
        print(f"  > Decay width: {width_nc:.6e} GeV")
    except (ValueError, TypeError) as ex:
        print(f"  > Calculation failed: {ex}")
