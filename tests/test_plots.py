# Importing the module containing the plotting functions to be tested
import pytest
import numpy as np
import matplotlib.pyplot as plt
from unittest import mock
from wt import postprocessing as pp  

# Fixture to prevent plt.show() from blocking tests
@pytest.fixture(autouse=True)
def prevent_show(monkeypatch):
    """Prevent plt.show() from blocking tests."""
    # Replace plt.show() with a no-op to avoid GUI pop-ups during testing
    monkeypatch.setattr(plt, "show", lambda: None)

# Test for the airfoil shape plotting function
def test_plot_airfoil_shape():
    # Generate x and y coordinates for a simple sinusoidal airfoil shape
    x = np.linspace(0, 1, 50)
    y = np.sin(2 * np.pi * x) * 0.1
    # Call the function to plot the airfoil shape
    pp.plot_airfoil_shape(x, y, airfoil_id="TestAirfoil", color="green", linestyle="--")
    # If no exception is raised, the test passes

# Test for plotting multiple airfoils
def test_plot_all_airfoils():
    # Create a dictionary containing data for two airfoils
    airfoil_data = {
        "001": {"x_coords": np.linspace(0, 1, 50), "y_coords": np.sin(2 * np.pi * np.linspace(0, 1, 50)) * 0.1},
        "002": {"x_coords": np.linspace(0, 1, 50), "y_coords": np.cos(2 * np.pi * np.linspace(0, 1, 50)) * 0.1}
    }
    # Call the function to plot all airfoils
    pp.plot_all_airfoils(airfoil_data)
    # If no exception is raised, the test passes

# Test for plotting power and thrust curves
def test_plot_power_thrust_curves():
    # Create a dictionary containing operational data for a wind turbine
    operational_data = {
        # Wind speeds in m/s
        "wind_speed_ms": np.array([3, 5, 8, 12, 15, 20]), 
        # Pitch angles in degrees 
        "pitch_deg": np.array([2, 1.5, 1, 0.5, 0.2, 0]),
        # Rotor speeds in RPM  
        "rot_speed_rpm": np.array([10, 20, 30, 40, 50, 60]),  
        # Aerodynamic power in kW
        "aero_power_kw": np.array([0, 100, 500, 1000, 1500, 2000]),  
        # Aerodynamic thrust in kW
        "aero_thrust_kw": np.array([1, 2, 3, 4, 5, 6])  
    }
    # Call the function to plot power and thrust curves
    pp.plot_power_thrust_curves(operational_data)
    # If no exception is raised, the test passes
