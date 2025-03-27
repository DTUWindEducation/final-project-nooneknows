import numpy as np
import os
import pandas as pd 
def load_blade_geometry(filepath="./inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat"):
    """Load blade geometry from AeroDyn15 blade definition file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
        index = next(i for i, line in enumerate(lines) if "BlSpn" in line)

        header_line = lines[index].strip()  # Get the second line (index 1)
        names = header_line.split()  # Split by whitespace

        units = lines[index+1].strip()
        unit = units.split()

    # Find where numerical data starts
    start_index = next(i for i, line in enumerate(lines) if "(m)" in line) 
    data = np.loadtxt(filepath, skiprows=start_index+1)

    # Create DataFrame for the data
    opt_data = pd.DataFrame(data, columns=names)
    # Optionally: Create a dictionary for the units
    units_opt = {names[i]: unit[i] for i in range(len(names))}


    



    # Read numerical data


    # Extract relevant columns
    r = data[:, 0]      # Spanwise position (m)
    beta = data[:, 4]   # Twist angle (degrees)
    chord = data[:, 5]  # Chord length (m)
    af_id = data[:, 6].astype(int)  # Airfoil ID (integer)

    return opt_data, units_opt

def load_airfoil_data(directory="./inputs/IEA-15-240-RWT/Airfoils"):
    """Load airfoil polars and coordinates from the airfoil directory."""
    airfoil_data = {}

    for af_id in range(0, 50):  # Airfoil IDs from 00 to 49
        af_key = f"{af_id:02d}"  # Format as "00", "01", ..., "49"

        # Load polar data
        polar_file = os.path.join(directory, f"IEA-15-240-RWT_AeroDyn15_Polar_{af_key}.dat")
        with open(polar_file, 'r') as f:
            lines = f.readlines()

        # Find where data starts
        start_index = next(i for i, line in enumerate(lines) if "!    (deg)" in line.lower())+2
        polar_data = np.loadtxt(polar_file, skiprows=start_index)

        alpha, Cl, Cd = polar_data[:, 0], polar_data[:, 1], polar_data[:, 2]

        # Load coordinates for airfoil shape
        coords_file = os.path.join(directory, f"IEA-15-240-RWT_AF{af_key}_Coords.txt")
        x_coords, y_coords = load_airfoil_coordinates(coords_file)

        # Store in dictionary
        airfoil_data[af_key] = {'alpha_deg': alpha, 'Cl': Cl, 'Cd': Cd, 'x_coords': x_coords, 'y_coords': y_coords}

    return airfoil_data

def load_airfoil_coordinates(file_path):
    """Load the airfoil coordinates (x/c, y/c) from the given file."""
    x_coords = []
    y_coords = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        start_index = next(i for i, line in enumerate(lines) if "!  x/c " in line) + 2

        for line in lines[start_index:]:
            if line.strip():
                parts = line.split()
                if len(parts) == 2:
                    x_coords.append(float(parts[0]))
                    y_coords.append(float(parts[1]))
                    
    return x_coords, y_coords

def load_operational_data(filepath="./inputs/IEA-15-240-RWT/IEA_15MW_RWT_Onshore.opt"):
    """Load operational data, such as wind speed, omega, and theta_p."""  
    # operational_data = {'V0': [], 'omega': {}, 'theta_p': {}}

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find where data starts
    start_index = next(i for i, line in enumerate(lines) if "wind speed [m/s]" in line.lower())+1
    operational_data = np.loadtxt(filepath, skiprows=start_index)

    wind_speed, pitch_deg, rot_speed_rpm, aero_power, aero_thrust = operational_data[:, 0], operational_data[:, 1], operational_data[:, 2] ,operational_data[:,3], operational_data[:,4]

    

    # Store in dictionary
    opt_data = {'wind_speed_ms': wind_speed, 'pitch_deg': pitch_deg, 'rot_speed_rpm': rot_speed_rpm, 'aero_power_kw': aero_power, 'aero_thrust_kw': aero_thrust}

                        
    return opt_data