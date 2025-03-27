import sys
import os

# Add the src directory to the Python path
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

# from wind_turbine_bem.bem_model import WindTurbineBEM
from load_files import load_blade_geometry, load_airfoil_data, load_operational_data

# Load data
blade_geometry = load_blade_geometry()
airfoil_data = load_airfoil_data()
operational_data = load_operational_data()

# Create and run BEM model
# bem_model = WindTurbineBEM(blade_geometry, airfoil_data, operational_data)
# bem_model.compute_power_curve