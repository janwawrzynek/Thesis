import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
csv_path = "rel_degree_of_Freedom.csv"
class RelDegreesOfFreedom:
    """
    A class to calculate the effective number of relativistic degrees of freedom, g_star,
    as a function of the dimensionless parameter z = m_N / T.
    """
    def __init__(self, csv_path, m_N):
        """
        Initializes the interpolator by converting temperature T from the data file
        into the dimensionless variable z.

        Args:
            csv_path (str): The file path to the CSV data for g_star vs. T.
            m_N (float): The mass of the particle N (e.g., a sterile neutrino) in GeV.
        """
        self.m_N = m_N
        df = pd.read_csv(csv_path)

        # Convert the temperature T from the file to z = m_N / T
        T_values = df["T (GeV)"]
        z_values = self.m_N / T_values

        # Get the g_star values
        g_star_values = df["g∗s"]

        # Since T is decreasing, z = m_N / T is increasing.
        # No array reversal is needed for interp1d.
        logZ = np.log10(z_values.to_numpy())
        logg = np.log10(g_star_values.to_numpy())

        # Create the interpolation function based on log10(z)
        self._interp = interp1d(logZ, logg, kind="linear", fill_value="extrapolate")

    def g_star(self, z):
        """
        Calculates g_star for a given z.

        Args:
            z (float or np.ndarray): The dimensionless value(s) z = m_N / T.

        Returns:
            float or np.ndarray: The value(s) of g_star.
        """
        # The interpolation is directly in terms of log10(z)
        return 10 ** self._interp(np.log10(z))


