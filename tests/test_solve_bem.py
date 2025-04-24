
from pathlib import Path
import sys
import os
import numpy as np
from wt.solve_bem import BEMSolver

#Data are given to test the BEMSolver class 
blade_data = {
    'BlSpn': np.array([10, 20, 30]),  # Blade span
    'BlAFID': np.array([1, 1, 1]),    # Airfoil index
    'BlChord': np.array([2, 2.5, 3]), # Chord lenght
    'BlTwist': np.array([5, 7, 10])   # Twist angle
}

airfoil_data = {
    '00': {
        'alpha_deg': np.array([-5, 0, 5, 10, 15]), #Angle of attack
        'Cl': np.array([0.5, 1.0, 1.2, 1.3, 1.1]), #Lift coefficient
        'Cd': np.array([0.01, 0.02, 0.03, 0.05, 0.07]) #Drag coefficient
    }
}

operational_data = {
    'wind_speed_ms': np.array([5, 10, 15]), #Wind speed
    'pitch_deg': np.array([0, 2, 5]), #Blade pitch angle
    'rot_speed_rpm': np.array([10, 15, 20]),  #Rotational speed
    'aero_power_kw': np.array([1000, 5000, 15678]) #Aerodynamic power
}

#From the BEMsolver class
bem = BEMSolver(blade_data, airfoil_data, operational_data)

def test_compute_tsr():
    #given
    wind_speed = 10  
    rot_speed_rpm = 15  
    
    # Calculate the expected tip speed ratio
    expected_tsr = (rot_speed_rpm * 2 * np.pi) / 60 * bem.R / wind_speed
    #when
    tsr = bem.compute_tsr(wind_speed, rot_speed_rpm)
    #then
    assert np.isclose(tsr, expected_tsr, rtol=1e-3)


def test_compute_aero_coeff():
    #given
    r = 20  #Span position
    alpha = 5  #Angle of attack
    
    # Calculate the coefficient of lift and drag
    Cl_exp, Cd_exp = bem.compute_aero_coeff(r, alpha)
    
    # Take the given data
    airfoil_id = '00'
    alpha_data = airfoil_data[airfoil_id]['alpha_deg']
    Cl_data = airfoil_data[airfoil_id]['Cl']
    Cd_data = airfoil_data[airfoil_id]['Cd']
    
    # Interpolate the coefficients
    Cl_int = np.interp(alpha, alpha_data, Cl_data)
    Cd_int = np.interp(alpha, alpha_data, Cd_data)
    
    assert np.isclose(Cl_exp, Cl_int, rtol=1e-3)
    assert np.isclose(Cd_exp, Cd_int, rtol=1e-3)


def test_get_optimal_operational_values():
    #given
    wind_speed = 10  # 
    
    # Calculate the expected optimal pitch and rotational speed
    optimal_pitch_exp = np.interp(wind_speed, operational_data['wind_speed_ms'], operational_data['pitch_deg'])
    optimal_rot_speed_exp = np.interp(wind_speed, operational_data['wind_speed_ms'], operational_data['rot_speed_rpm'])
    
    #when
    optimal_pitch, optimal_rot_speed = bem.get_optimal_operational_values(wind_speed)
    
    #then
    assert np.isclose(optimal_pitch, optimal_pitch_exp, rtol=1e-3)
    assert np.isclose(optimal_rot_speed, optimal_rot_speed_exp, rtol=1e-3)


def test_get_rated_wind_speed():
    
    # Calculate the expected rated wind speed
    rated_ws_exp = operational_data['wind_speed_ms'][np.argmax(operational_data['rot_speed_rpm'])]
    
    # when
    rated_ws = bem.get_rated_wind_speed()
    
    # then
    assert np.isclose(rated_ws, rated_ws_exp, rtol=1e-3)
