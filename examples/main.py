import sys
import os
import matplotlib.pyplot as plt

# Add the src directory to the Python path
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

# from wind_turbine_bem.bem_model import WindTurbineBEM
from load_files import load_blade_geometry, load_airfoil_data, load_operational_data
from postprocessing import plot_all_airfoils

# Load data
opt_data, units_opt = load_blade_geometry()
airfoil_data = load_airfoil_data()
operational_data = load_operational_data()
plot_all_airfoils(airfoil_data)

# Create and run BEM model
# bem_model = WindTurbineBEM(blade_geometry, airfoil_data, operational_data)
# bem_model.compute_power_curve