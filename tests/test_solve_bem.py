"""
Unit tests for the BEMSolver class.

These tests validate the functionality of the BEMSolver methods, including
tip speed ratio computation, aerodynamic coefficient interpolation, and
optimal operational value determination.
"""

import numpy as np
from wt.solve_bem import BEMSolver

# Data for testing the BEMSolver class
blade_data = {
    'BlSpn': np.array([10, 20, 30]),  # Blade span
    'BlAFID': np.array([1, 1, 1]),    # Airfoil index
    'BlChord': np.array([2, 2.5, 3]),  # Chord length
    'BlTwist': np.array([5, 7, 10])   # Twist angle
}

airfoil_data = {
    '00': {
        'alpha_deg': np.array([-5, 0, 5, 10, 15]),  # Angle of attack
        'Cl': np.array([0.5, 1.0, 1.2, 1.3, 1.1]),  # Lift coefficient
        'Cd': np.array([0.01, 0.02, 0.03, 0.05, 0.07])  # Drag coefficient
    }
}

operational_data = {
    'wind_speed_ms': np.array([5, 10, 15]),  # Wind speed
    'pitch_deg': np.array([0, 2, 5]),  # Blade pitch angle
    'rot_speed_rpm': np.array([10, 15, 20]),  # Rotational speed
    'aero_power_kw': np.array([1000, 5000, 15678])  # Aerodynamic power
}

# Initialize the BEMSolver class
bem = BEMSolver(blade_data, airfoil_data, operational_data)


def test_compute_tsr():
    """
    Test the compute_tsr method for calculating the tip speed ratio.
    """
    wind_speed = 10
    rot_speed_rpm = 15

    # Calculate the expected tip speed ratio
    expected_tsr = (rot_speed_rpm * 2 * np.pi) / 60 * bem.R / wind_speed

    # Compute the tip speed ratio
    tsr = bem.compute_tsr(wind_speed, rot_speed_rpm)

    # Assert the result is close to the expected value
    assert np.isclose(tsr, expected_tsr, rtol=1e-3)


def test_compute_aero_coeff():
    """
    Test the compute_aero_coeff method for interpolating aerodynamic
    coefficients.
    """
    r = 20  # Span position
    alpha = 5  # Angle of attack

    # Compute the aerodynamic coefficients
    Cl_exp, Cd_exp = bem.compute_aero_coeff(r, alpha)

    # Extract airfoil data
    airfoil_id = '00'
    alpha_data = airfoil_data[airfoil_id]['alpha_deg']
    Cl_data = airfoil_data[airfoil_id]['Cl']
    Cd_data = airfoil_data[airfoil_id]['Cd']

    # Interpolate the coefficients
    Cl_int = np.interp(alpha, alpha_data, Cl_data)
    Cd_int = np.interp(alpha, alpha_data, Cd_data)

    # Assert the results are close to the interpolated values
    assert np.isclose(Cl_exp, Cl_int, rtol=1e-3)
    assert np.isclose(Cd_exp, Cd_int, rtol=1e-3)


def test_get_optimal_operational_values():
    """
    Test the get_optimal_operational_values method for determining the
    optimal pitch and rotational speed for a given wind speed.
    """
    wind_speed = 10

    # Calculate the expected optimal pitch and rotational speed
    optimal_pitch_exp = np.interp(
        wind_speed, operational_data['wind_speed_ms'],
        operational_data['pitch_deg']
    )
    optimal_rot_speed_exp = np.interp(
        wind_speed, operational_data['wind_speed_ms'],
        operational_data['rot_speed_rpm']
    )

    # Compute the optimal values
    optimal_pitch, optimal_rot_speed = bem.get_optimal_operational_values(
        wind_speed
    )

    # Assert the results are close to the expected values
    assert np.isclose(optimal_pitch, optimal_pitch_exp, rtol=1e-3)
    assert np.isclose(optimal_rot_speed, optimal_rot_speed_exp, rtol=1e-3)


def test_get_rated_wind_speed():
    """
    Test the get_rated_wind_speed method for determining the rated wind speed.
    """
    # Calculate the expected rated wind speed
    rated_ws_exp = operational_data['wind_speed_ms'][
        np.argmax(operational_data['rot_speed_rpm'])
    ]

    # Compute the rated wind speed
    rated_ws = bem.get_rated_wind_speed()

    # Assert the result is close to the expected value
    assert np.isclose(rated_ws, rated_ws_exp, rtol=1e-3)
