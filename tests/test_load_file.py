"""
Unit tests for the load_files module.

These tests validate the functionality of the loaders for operational data,
blade geometry, and airfoil data.
"""

from pathlib import Path
import numpy as np
from wt import load_files

# Define the path to the inputs directory relative to the test file
DATA_DIR = Path('./inputs/IEA-15-240-RWT')


def test_load_operational_data_default():
    """
    Test that OperationalDataLoader works as expected with default values.
    """
    path_load_operational_data_file = DATA_DIR / 'IEA_15MW_RWT_Onshore.opt'

    wind_speed_exp = 10
    aero_thrust_exp = 1898.680029

    operational_loader = load_files.OperationalDataLoader(
        path_load_operational_data_file
    )
    operational_data = operational_loader.load()

    wind_speed = operational_data['wind_speed_ms']
    aero_thrust = operational_data['aero_thrust_kw']

    assert np.isclose(wind_speed[6], wind_speed_exp, rtol=1e-3), (
        f"Expected wind speed {wind_speed_exp}, but got {wind_speed[6]}"
    )
    assert np.isclose(aero_thrust[6], aero_thrust_exp, rtol=1e-3), (
        f"Expected aero thrust {aero_thrust_exp}, but got {aero_thrust[6]}"
    )


def test_load_blade_geometry_default():
    """
    Test that BladeGeometryLoader works as expected with default values.
    """
    path_load_blade_geometry_file = (
        DATA_DIR / 'IEA-15-240-RWT_AeroDyn15_blade.dat'
    )

    bl_spn_exp = 4.775507409073585e+00

    blade_loader = load_files.BladeGeometryLoader(
        path_load_blade_geometry_file)
    opt_data, _ = blade_loader.load()

    bl_spn = opt_data['BlSpn']

    assert np.isclose(bl_spn[2], bl_spn_exp, rtol=1e-3), (
        f"Expected BlSpn {bl_spn_exp}, but got {bl_spn[2]}"
    )


def test_load_airfoil_data_default():
    """
    Test that AirfoilDataLoader works as expected with default values.
    """
    path_to_airfoil_dir = DATA_DIR / 'Airfoils'

    x_c_exp = 8.80309441682373e-01
    y_c_exp = 5.94483425110475e-03

    airfoil_loader = load_files.AirfoilDataLoader(
        directory=path_to_airfoil_dir)
    airfoil_data = airfoil_loader.load()

    af40_data = airfoil_data['40']

    x_c = af40_data['x_coords']
    y_c = af40_data['y_coords']

    assert np.isclose(x_c[5], x_c_exp, rtol=1e-3), (
        f"Expected x/c {x_c_exp}, but got {x_c[5]}"
    )
    assert np.isclose(y_c[5], y_c_exp, rtol=1e-3), (
        f"Expected y/c {y_c_exp}, but got {y_c[5]}"
    )
