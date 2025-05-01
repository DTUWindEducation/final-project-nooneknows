"""
Main script for performing Blade Element Momentum (BEM) calculations.

This script loads blade geometry, airfoil, and operational data, initializes
the BEM solver, and performs calculations for thrust, torque, power, and
aerodynamic coefficients. It also generates plots for visualization.
"""

import numpy as np
from pathlib import Path
from wt import load_files, solve_bem, postprocessing

# Get the directory of this script regardless of where it is run
THIS_DIR = Path(__file__).resolve().parent

# Define paths relative to this script's location
BASE_PATH = THIS_DIR / ".." / "inputs" / "IEA-15-240-RWT"
BLADE_FILE = BASE_PATH / "IEA-15-240-RWT_AeroDyn15_blade.dat"
AIRFOIL_DIR = BASE_PATH / "Airfoils"
OPERATIONAL_FILE = BASE_PATH / "IEA_15MW_RWT_Onshore.opt"

# Load blade geometry data
blade_loader = load_files.BladeGeometryLoader(filepath=BLADE_FILE)
blade_data, _ = blade_loader.load()

# Load airfoil data
airfoil_loader = load_files.AirfoilDataLoader(directory=AIRFOIL_DIR)
airfoil_data = airfoil_loader.load()

# Load operational data
operational_loader = load_files.OperationalDataLoader(filepath=OPERATIONAL_FILE)
operational_data = operational_loader.load()

# Initialize the BEM solver with the loaded data
bem_solver = solve_bem.BEMSolver(blade_data, airfoil_data, operational_data)

# Perform an initial BEM calculation for specific conditions
# Example values: wind speed = 10 m/s, pitch angle = 0 degrees, rotor speed = 5 rpm
thrust, torque, power, CT, CP = bem_solver.solve_bem(10, 0, 5)

# Plot all airfoil data for visualization
postprocessing.plot_all_airfoils(airfoil_data)

# Define specific conditions for BEM calculation
WIND_SPEED = 10  # Wind speed in m/s
PITCH_ANGLE = 0  # Pitch angle in degrees
ROT_SPEED_RPM = 5  # Rotor speed in revolutions per minute (rpm)

# Perform BEM calculation for the defined conditions
thrust, torque, power, CT, CP = bem_solver.solve_bem(
    WIND_SPEED, PITCH_ANGLE, ROT_SPEED_RPM
)

# Print the results of the BEM calculation
print(f"Thrust: {thrust:.2f} N")  # Thrust force in Newtons
print(f"Torque: {torque:.2f} Nm")  # Torque in Newton-meters
print(f"Power: {power:.2f} W")  # Power in Watts
print(f"Thrust Coefficient (CT): {CT:.4f}")  # Non-dimensional thrust coefficient
print(f"Power Coefficient (CP): {CP:.4f}")  # Non-dimensional power coefficient

# Plot power and thrust curves based on operational data
postprocessing.plot_power_thrust_curves(operational_data)

# Generate a range of wind speeds for power and thrust curve calculations
WIND_SPEEDS = np.linspace(3, 25, 100)  # Wind speeds from 3 m/s to 25 m/s

# Compute power and thrust curves for the range of wind speeds
power_curve, thrust_curve, omega, pitch = bem_solver.compute_power_thrust_curve(WIND_SPEEDS)

# Plot the computed power and thrust curves
postprocessing.plot_bem_power_thrust_curves(WIND_SPEEDS, power_curve, thrust_curve)
postprocessing.plot_pitch_rot_speed(WIND_SPEEDS, pitch, omega)

wind_speed = 10.0
pitch, rpm, _ = bem_solver.get_optimal_operational_values(wind_speed)

# Compute and plot spanwise normal and tangential loads
spanwise_data = bem_solver.compute_spanwise_normal_tangential_loads(wind_speed, pitch, rpm)
postprocessing.plot_spanwise_normal_tangential_loads(spanwise_data)







