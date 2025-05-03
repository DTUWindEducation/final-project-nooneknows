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

""" 
    1. Load and parse the provided turbine data
"""

# Load blade geometry data
blade_loader = load_files.BladeGeometryLoader(filepath=BLADE_FILE)
blade_data, _ = blade_loader.load()

# Load airfoil data
airfoil_loader = load_files.AirfoilDataLoader(directory=AIRFOIL_DIR)
airfoil_data = airfoil_loader.load()

# Load operational data
operational_loader = load_files.OperationalDataLoader(
    filepath=OPERATIONAL_FILE)
operational_data = operational_loader.load()

"""
    3.Compute lift coefficient (Cl) and drag coefficient (Cd) as function of span position (r) and angle of attack (α)
"""
# Choose a representative span position (e.g., mid-span) and angle of attack
r_span = blade_data['BlSpn'][len(blade_data['BlSpn']) // 2]  # Mid-span
alpha_deg = 5.0  # Example angle of attack in degrees
bem_solver = solve_bem.BEMSolver(blade_data, airfoil_data, operational_data)
# Compute lift and drag coefficients at that position
Cl, Cd = bem_solver.compute_aero_coeff(r_span, alpha_deg)
print(f"At r = {r_span:.2f} m and alpha = {alpha_deg}° --> Cl = {Cl:.4f}, Cd = {Cd:.4f}")

"""
    4.Compute the axial (a) and tangential (a′) induction factors as function of span position (r),
    the inflow wind speed V0, the blade pitch angle (θp) and the rotational speed ω
"""
# Define parameters
ROT_SPEED_RPM  = 7.101976
omega_rad_s = (ROT_SPEED_RPM * 2 * np.pi) / 60
chord = blade_data['BlChord'][len(blade_data['BlChord']) // 2]
twist_deg = blade_data['BlTwist'][len(blade_data['BlTwist']) // 2]
twist_rad = np.radians(twist_deg)
airfoil_id = str(int(blade_data['BlAFID'][len(blade_data['BlAFID']) // 2]) - 1).zfill(2)
WIND_SPEED = 10.000000 # [m/s]
PITCH_ANGLE = 0.000535 # [deg]

# Compute induction factors
a, a_prime = bem_solver.compute_induction(
    r_span, WIND_SPEED, PITCH_ANGLE, omega_rad_s, chord, twist_rad, airfoil_id
)
print(f"At r = {r_span:.2f} m --> a = {a:.4f}, a' = {a_prime:.4f}")
"""
    5.Compute the thrust (T), torque (M), and power (P) of the rotor as function of the inflow wind speed V0 ,
    the blade pitch angle (θp) and the rotational speed ω
"""

# Initialize the BEM solver with the loaded data
bem_solver = solve_bem.BEMSolver(blade_data, airfoil_data, operational_data)

# Perform an initial BEM calculation for specific conditions
# Example values: wind speed = 10 m/s,
# pitch angle = 0 degrees, rotor speed = 5 rpm
thrust, torque, power, CT, CP = bem_solver.solve_bem(10, 0, 5)

"""
    2. Plot the provided airfoil shapes in one figure
"""

# Plot all airfoil data for visualization
postprocessing.plot_all_airfoils(airfoil_data)

"""
    6. Compute optimal operational strategy, i.e., blade pitch angle (θp) and rotational speed (ω),
    as function of wind speed (V0), based on the provided operational strategy in IEA_15MW_RWT_Onshore.opt
"""
# Define specific conditions for BEM calculation
WIND_SPEED = 10  # Wind speed in m/s
PITCH_ANGLE = 0  # Pitch angle in degrees
ROT_SPEED_RPM = 5  # Rotor speed in revolutions per minute (rpm)

# Perform BEM calculation for the defined conditions
thrust, torque, power, CT, CP = bem_solver.solve_bem(
    WIND_SPEED, PITCH_ANGLE, ROT_SPEED_RPM
)

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

"""
    7. Compute and plot power curve ($P(V_0)$) and thrust curve ($T(V_0)$) 
    based on the optimal operational strategy obtained in the previous function.
"""

# Plot power and thrust curves from operational data file
postprocessing.plot_power_thrust_curves(operational_data)

# Generate a range of wind speeds for power and thrust curve calculations
WIND_SPEEDS = np.linspace(3, 25, 100)  # Wind speeds from 3 m/s to 25 m/s

# Compute power and thrust curves for the range of wind speeds
power_curve, thrust_curve, omega, pitch = bem_solver.compute_power_thrust_curve(WIND_SPEEDS)

# Plot the computed power and thrust curves
postprocessing.plot_bem_power_thrust_curves(
    WIND_SPEEDS, power_curve, thrust_curve)

"""
EXTRA FUNCTION n.1: Plot the optimal pitch angle and rotational speed curves
"""

postprocessing.plot_pitch_rot_speed(WIND_SPEEDS, pitch, omega)

wind_speed = 10.0
pitch, rpm, _ = bem_solver.get_optimal_operational_values(wind_speed)

"""
EXTRA FUNCTION n.2: Compute spanwise normal and tangential loads 
EXTRA FUNCTION n.3: Plot spanwise normal and tangential loads
"""

# Compute and plot spanwise normal and tangential loads
spanwise_data = bem_solver.compute_spanwise_normal_tangential_loads(
    wind_speed, pitch, rpm)
postprocessing.plot_spanwise_normal_tangential_loads(spanwise_data)
