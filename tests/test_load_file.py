

# Import necessary modules
from pathlib import Path
import numpy as np
from wt import load_files

# Define the path to the inputs directory relative to the test file
DATA_DIR = Path('./inputs/IEA-15-240-RWT')

def test_load_OperationalData_default():
    """Test that load_files works as expected with default values."""

    # Define the path to the operational data file
    path_load_OperationalData_file = DATA_DIR / 'IEA_15MW_RWT_Onshore.opt'

    # Define the expected wind speed and aero thrust values at index 6
    wind_speed_exp = 10
    aero_thrust_exp = 1898.680029

    # Create an instance of OperationalDataLoader and load the data
    operational_loader = load_files.OperationalDataLoader(path_load_OperationalData_file)
    operational_data = operational_loader.load()

    # Extract wind speed and aero thrust data from the loaded dictionary
    wind_speed = operational_data['wind_speed_ms']
    aero_thrust = operational_data['aero_thrust_kw']

    # Assert that the wind speed matches the expected value
    assert np.isclose(wind_speed[6], wind_speed_exp, rtol=1e-3)

    # Assert that the aero thrust matches the expected value
    assert np.isclose(aero_thrust[6], aero_thrust_exp, rtol=1e-3)

def test_load_BladeGeometry_default():
    """Test that load_files works as expected with default values."""

    # Define the path to the blade geometry file
    path_load_BladeGeometry_file = DATA_DIR / 'IEA-15-240-RWT_AeroDyn15_blade.dat'

    # Define the expected BlSpn value at index 3
    BlSpn_exp = 4.775507409073585e+00

    # Create an instance of BladeGeometryLoader and load the data
    blade_loader = load_files.BladeGeometryLoader(path_load_BladeGeometry_file)
    opt_data, _ = blade_loader.load()

    # Extract the BlSpn column from the loaded data
    BlSpn = opt_data['BlSpn']

    # Assert that the BlSpn value matches the expected value
    assert np.isclose(BlSpn[2], BlSpn_exp, rtol=1e-3)

def test_load_AirfoilData_default():
    """Test that AirfoilDataLoader works as expected with default values."""

    # Define the path to the Airfoils directory
    path_to_airfoil_dir = DATA_DIR / 'Airfoils'

    # Define the expected airfoil coordinates at index 5
    x_c_exp = 8.80309441682373e-01
    y_c_exp = 5.94483425110475e-03

    # Create an instance of AirfoilDataLoader and load the data
    airfoil_loader = load_files.AirfoilDataLoader(directory=path_to_airfoil_dir)
    airfoil_data = airfoil_loader.load()

    # Access the data for airfoil '40' (AF40)
    af40_data = airfoil_data['40']

    # Extract x/c and y/c coordinates from the airfoil data
    x_c = af40_data['x_coords']
    y_c = af40_data['y_coords']

    # Assert that the x/c value matches the expected value
    assert np.isclose(x_c[5], x_c_exp, rtol=1e-3), f"Expected x/c {x_c_exp}, but got {x_c[5]}"

    # Assert that the y/c value matches the expected value
    assert np.isclose(y_c[5], y_c_exp, rtol=1e-3), f"Expected y/c {y_c_exp}, but got {y_c[5]}"

