# test_load_file.py
from pathlib import Path
import sys
import os
import numpy as np

# Add the path to the src folder to import load_resp from src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

# Import the function load_files from load_resp module
import load_files

# Define the path to your inputs directory (relative to the test file)
DATA_DIR = Path(__file__).resolve().parent.parent / 'inputs' / 'IEA-15-240-RWT'

def test_load_OperationalData_default():
    """Test that load_files works as expected with default values."""
    
    # Path to the input file to load
    path_load_OperationalData_file = DATA_DIR / 'IEA_15MW_RWT_Onshore.opt'
    
    # Expected values for wind speed and aero thrust at index 6
    wind_speed_exp = 10  # Expected wind speed value
    aero_thrust_exp = 1898.680029  # Expected aero thrust value
    
    # Create the OperationalDataLoader object and load the data
    operational_loader = load_files.OperationalDataLoader(path_load_OperationalData_file)
    operational_data = operational_loader.load()
    
    # Extract the relevant data from the dictionary
    wind_speed = operational_data['wind_speed_ms']
    aero_thrust = operational_data['aero_thrust_kw']
    
    # Assert that the wind speed and aero thrust values match the expected values
    assert np.isclose(wind_speed[6], wind_speed_exp, rtol=1e-3)
    
    assert np.isclose(aero_thrust[6], aero_thrust_exp, rtol=1e-3)


def test_load_BladeGeometry_default():
    """Test that load_files works as expected with default values."""
    
    # Path to the input file to load
    path_load_BladeGeometry_file = DATA_DIR / 'IEA-15-240-RWT_AeroDyn15_blade.dat'
    
    # Expected values for BlSpn at index 3
    BlSpn_exp  = 4.775507409073585e+00  # Expected value
     
    # Create the BladeGeometryLoader object and load the data
    blade_loader = load_files.BladeGeometryLoader(path_load_BladeGeometry_file)
    opt_data, _ = blade_loader.load()
    
    # Extract the BlSpn column
    BlSpn = opt_data['BlSpn']
    
    # Assert that the BlSpn value matches the expected value
    assert np.isclose(BlSpn[2], BlSpn_exp, rtol=1e-3)
    
def test_load_AirfoilData_default():
    """Test that AirfoilDataLoader works as expected with default values."""

    # Path to the Airfoils directory
    path_to_airfoil_dir = DATA_DIR / 'Airfoils'
    
    # Expected values for airfoil coordinates at index 5
    x_c_exp  = 8.80309441682373e-01  # Expected value
    y_c_exp  = 5.94483425110475e-03  # Expected value
     
    # Create the AirfoilDataLoader object and load the data
    airfoil_loader = load_files.AirfoilDataLoader(directory=path_to_airfoil_dir)
    airfoil_data = airfoil_loader.load()
    
    # Access airfoil '40' (AF40)
    af40_data = airfoil_data['40']
    
    x_c = af40_data['x_coords']
    y_c = af40_data['y_coords']
    
    # Assert that the x/c and y/c values match expected values
    assert np.isclose(x_c[5], x_c_exp, rtol=1e-3), f"Expected x/c {x_c_exp}, but got {x_c[5]}"
    assert np.isclose(y_c[5], y_c_exp, rtol=1e-3), f"Expected y/c {y_c_exp}, but got {y_c[5]}"

