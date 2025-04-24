# Import necessary libraries and modules
import numpy as np
from pathlib import Path
from wt import load_files, solve_bem, postprocessing

# Get the directory of this script regardless of where it is run
this_dir = Path(__file__).resolve().parent

# Define paths relative to this script's location

# These paths point to input files required for the simulation
#general folder
base_path = this_dir / ".." / "inputs" / "IEA-15-240-RWT" 
# Blade geometry data
blade_file = base_path / "IEA-15-240-RWT_AeroDyn15_blade.dat" 
# Directory containing airfoil data
airfoil_dir = base_path / "Airfoils"  
# Operational data file
operational_file = base_path / "IEA_15MW_RWT_Onshore.opt"  

# Initialize loaders with pathlib paths
# Load blade geometry data
blade_loader = load_files.BladeGeometryLoader(filepath=blade_file)
blade_data, _ = blade_loader.load()

# Load airfoil data
airfoil_loader = load_files.AirfoilDataLoader(directory=airfoil_dir)
airfoil_data = airfoil_loader.load()

# Load operational data
operational_loader = load_files.OperationalDataLoader(filepath=operational_file)
operational_data = operational_loader.load()

# Initialize the BEM (Blade Element Momentum) solver with the loaded data
bem_solver = solve_bem.BEMSolver(blade_data, airfoil_data, operational_data)

# Perform an initial BEM calculation for a specific wind speed, pitch angle, and rotor speed
# Example values: wind speed = 10 m/s, pitch angle = 0 degrees, rotor speed = 5 rpm
thrust, torque, power, CT, CP = bem_solver.solve_bem(10, 0, 5)

# Plot all airfoil data for visualization
postprocessing.plot_all_airfoils(airfoil_data)

# Define specific conditions for BEM calculation
# Wind speed in m/s
wind_speed = 10 
# Pitch angle in degrees
pitch_angle = 0  
# Rotor speed in revolutions per minute (rpm)
rot_speed_rpm = 5  

# Perform BEM calculation for the defined conditions
thrust, torque, power, CT, CP = bem_solver.solve_bem(wind_speed, pitch_angle, rot_speed_rpm)

# Print the results of the BEM calculation
# Thrust force in Newtons
print(f"Thrust: {thrust:.2f} N")  
# Torque in Newton-meters
print(f"Torque: {torque:.2f} Nm")  
# Power in Watts
print(f"Power: {power:.2f} W")  
# Non-dimensional thrust coefficient
print(f"Thrust Coefficient (CT): {CT:.4f}")  
 # Non-dimensional power coefficient
print(f"Power Coefficient (CP): {CP:.4f}") 

# Plot power and thrust curves based on operational data
postprocessing.plot_power_thrust_curves(operational_data)

# Generate a range of wind speeds for power and thrust curve calculations
# Wind speeds from 3 m/s to 25 m/s
wind_speeds = np.linspace(3, 25, 100)  

# Compute power and thrust curves for the range of wind speeds
power_curve, thrust_curve = bem_solver.compute_power_thrust_curve(wind_speeds)

# Plot the computed power and thrust curves
bem_solver.plot_power_thrust_curves(wind_speeds, power_curve, thrust_curve)

