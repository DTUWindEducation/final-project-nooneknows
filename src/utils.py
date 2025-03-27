import numpy as np
import os

def load_blade_geometry(filepath="./inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat"):
    """Load blade geometry from AeroDyn15 blade definition file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find where numerical data starts
    start_index = next(i for i, line in enumerate(lines) if "BlSpn" in line) + 2

    # Read numerical data
    data = np.loadtxt(filepath, skiprows=start_index)

    # Extract relevant columns
    r = data[:, 0]      # Spanwise position (m)
    beta = data[:, 4]   # Twist angle (degrees)
    chord = data[:, 5]  # Chord length (m)
    af_id = data[:, 6].astype(int)  # Airfoil ID (integer)

    return {'r': r, 'chord': chord, 'beta': beta, 'af_id': [f"{af:02d}" for af in af_id]}

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
        airfoil_data[af_key] = {'alpha': alpha, 'Cl': Cl, 'Cd': Cd, 'x_coords': x_coords, 'y_coords': y_coords}

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
    operational_data = {'V0': [], 'omega': {}, 'theta_p': {}}
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Read the file using pandas, assuming it's space-separated
    data = pd.read_csv(filepath, delim_whitespace=True, header=None)

    # Define column names based on the provided data
    data.columns = ['wind_speed', 'pitch', 'rot_speed', 'aero_power', 'aero_thrust']

    # Display the first few rows of the data
    print(data)
                        


    return operational_data