"""
Module for Blade Element Momentum (BEM) calculations.

This module defines the BEMSolver class, which provides methods for computing
aerodynamic coefficients, induction factors, and power/thrust curves for wind
turbines.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class BEMSolver:
    """
    A class to handle Blade Element Momentum (BEM) calculations.
    """

    def __init__(self, blade_data, airfoil_data, operational_data, B=3, rho=1.225):
        """
        Initialize the BEMSolver with blade, airfoil, and operational data.

        Args:
            blade_data (dict): Blade geometry data.
            airfoil_data (dict): Airfoil aerodynamic data.
            operational_data (dict): Operational data.
            B (int): Number of blades (default: 3).
            rho (float): Air density in kg/m^3 (default: 1.225).
        """
        self.blade_data = blade_data
        self.airfoil_data = airfoil_data
        self.operational_data = operational_data
        self.B = B
        self.rho = rho
        self.R = blade_data['BlSpn'].max()
        self.A = np.pi * self.R**2

    def compute_tsr(self, wind_speed, rot_speed_rpm):
        """
        Compute the tip speed ratio (TSR).

        Args:
            wind_speed (float): Wind speed in m/s.
            rot_speed_rpm (float): Rotational speed in RPM.

        Returns:
            float: Tip speed ratio.
        """
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        return omega * self.R / wind_speed

    def compute_aero_coeff(self, r, alpha):
        """
        Compute lift and drag coefficients for a given span position and angle
        of attack.

        Args:
            r (float): Span position in meters.
            alpha (float): Angle of attack in degrees.

        Returns:
            tuple: Lift coefficient (Cl) and drag coefficient (Cd).
        """
        idx = (np.abs(self.blade_data['BlSpn'] - r)).argmin()
        airfoil_id = str(int(self.blade_data['BlAFID'][idx]) - 1).zfill(2)
        alpha_data = self.airfoil_data[airfoil_id]['alpha_deg']
        Cl_data = self.airfoil_data[airfoil_id]['Cl']
        Cd_data = self.airfoil_data[airfoil_id]['Cd']
        Cl = np.interp(alpha, alpha_data, Cl_data)
        Cd = np.interp(alpha, alpha_data, Cd_data)
        return Cl, Cd

    def compute_induction(self, r, V0, theta_p, omega, c_r, beta_r, airfoil_id):
        """
        Compute axial (a) and tangential (a') induction factors.

        Args:
            r (float): Span position in meters.
            V0 (float): Wind speed in m/s.
            theta_p (float): Pitch angle in degrees.
            omega (float): Rotational speed in rad/s.
            c_r (float): Chord length in meters.
            beta_r (float): Twist angle in radians.
            airfoil_id (str): Airfoil ID.

        Returns:
            tuple: Axial induction factor (a) and tangential induction factor (a').
        """
        a, a_prime = 0, 0
        solidity = (c_r * self.B) / (2 * np.pi * r + 1e-6)

        for _ in range(100):
            phi = np.arctan((1 - a) * V0 / ((1 + a_prime) * omega * r + 1e-6))
            alpha = np.degrees(phi) - (theta_p + np.rad2deg(beta_r))
            Cl, Cd = self.compute_aero_coeff(r, alpha)
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
            a_new = 1 / (4 * np.sin(phi)**2 / (solidity * Cn) + 1)
            a_prime_new = 1 / (4 * np.sin(phi) * np.cos(phi) /
                               (solidity * Ct) - 1)
            if np.abs(a - a_new) < 1e-4 and np.abs(a_prime - a_prime_new) < 1e-4:
                break
            a, a_prime = a_new, a_prime_new

        return a_new, a_prime_new

    def solve_bem(self, V0, theta_p, rot_speed_rpm, target_power_kw = None):
        """
        Solve BEM equations to compute thrust, torque, and power.

        Args:
            V0 (float): Wind speed in m/s.
            theta_p (float): Pitch angle in degrees.
            rot_speed_rpm (float): Rotational speed in RPM.

        Returns:
            tuple: Thrust (N), torque (Nm), power (W), thrust coefficient (CT),
                   and power coefficient (CP).
        """
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        r_vals = self.blade_data['BlSpn']
        thrust, torque = 0, 0
        # rated_power = 15678.725250 * 1e3
        dT_vals, dM_vals = [], []

        for i, r in enumerate(self.blade_data['BlSpn']):
            c_r = self.blade_data['BlChord'][i]
            beta_r = np.radians(self.blade_data['BlTwist'][i])
            airfoil_id = str(int(self.blade_data['BlAFID'][i]) - 1).zfill(2)
            a, a_prime = self.compute_induction(r, V0, theta_p, omega, c_r,
                                                beta_r, airfoil_id)
            dr = self.R / len(self.blade_data)
            dT = 4 * np.pi * r * self.rho * V0**2 * a * (1 - a) * dr
            dM = 4 * np.pi * r**3 * self.rho * V0 * omega * a_prime * (1 - a) * dr
            dT_vals.append(dT)
            dM_vals.append(dM)

        thrust = np.trapz(dT_vals, dx=dr)
        torque = np.trapz(dM_vals, dx=dr)
        power = torque * omega

        if target_power_kw is not None:
            target_power_kw=target_power_kw*1e3
            power = min(power,target_power_kw)
        
        CT = thrust / (0.5 * self.rho * self.A * V0**2)
        CP = power / (0.5 * self.rho * self.A * V0**3)

        return thrust, torque, power, CT, CP

    def get_optimal_operational_values(self, wind_speed):
        """
        Get optimal pitch angle and rotational speed for a given wind speed.

        Args:
            wind_speed (float): Wind speed in m/s.

        Returns:
            tuple: Optimal pitch angle (degrees) and rotational speed (RPM).
        """
        wind_speeds = self.operational_data['wind_speed_ms']
        optimal_pitch_angles = self.operational_data['pitch_deg']
        optimal_rot_speeds = self.operational_data['rot_speed_rpm']
        power_curve_kw = self.operational_data['aero_power_kw']

        # Interpolate pitch, rotation speed, and power
        interp_pitch = np.interp(wind_speed, wind_speeds, optimal_pitch_angles)
        interp_rot_speed = np.interp(wind_speed, wind_speeds, optimal_rot_speeds)
        interp_power = np.interp(wind_speed, wind_speeds, power_curve_kw)

        # Identify rated power and rated wind speed
        rated_power = np.max(power_curve_kw)
        rated_indices = np.where(power_curve_kw == rated_power)[0]
        rated_index = rated_indices[0]
        rated_wind_speed = wind_speeds[rated_index]
        rated_rot_speed = optimal_rot_speeds[rated_index]

        # If wind speed >= rated, cap at rated power behavior
        if wind_speed >= rated_wind_speed:
            optimal_pitch = interp_pitch   # Pitch still varies to regulate power
            optimal_rot_speed = rated_rot_speed  # Fix RPM at rated
        else:
            optimal_pitch = interp_pitch
            optimal_rot_speed = interp_rot_speed

        return optimal_pitch, optimal_rot_speed, power_curve_kw

    def get_rated_wind_speed(self):
        """
        Get the rated wind speed based on maximum rotational speed.

        Returns:
            float: Rated wind speed in m/s.
        """
        max_rot_speed = np.max(self.operational_data['rot_speed_rpm'])
        rated_index = np.where(self.operational_data['rot_speed_rpm'] ==
                               max_rot_speed)[0][0]
        return self.operational_data['wind_speed_ms'][rated_index]

    def compute_power_thrust_curve(self, wind_speeds):
        """
        Compute power and thrust curves for a range of wind speeds.

        Args:
            wind_speeds (np.ndarray): Array of wind speeds in m/s.

        Returns:
            tuple: Power curve (W) and thrust curve (N).
        """
        power_curve = []
        thrust_curve = []
        omega_curve = []
        pitch_curve = []
        for V0 in wind_speeds:
            theta_p, omega_rpm, power_curve_kw = self.get_optimal_operational_values(V0)
            interp_power_kw = np.interp(V0, self.operational_data['wind_speed_ms'], power_curve_kw)
            T, M, P, CT, CP = self.solve_bem(V0, theta_p, omega_rpm, target_power_kw=interp_power_kw)
            power_curve.append(P)
            thrust_curve.append(T)
            omega_curve.append(omega_rpm)
            pitch_curve.append(theta_p)
        return np.array(power_curve), np.array(thrust_curve), np.array(omega_curve), np.array(pitch_curve)

    def plot_power_thrust_curves(self, wind_speeds, power_curve, thrust_curve):
        """
        Plot power and thrust curves as a function of wind speed.

        Args:
            wind_speeds (np.ndarray): Array of wind speeds in m/s.
            power_curve (np.ndarray): Power curve in W.
            thrust_curve (np.ndarray): Thrust curve in N.
        """
        fig, ax1 = plt.subplots()

        ax1.set_xlabel('Wind Speed (m/s)')
        ax1.set_ylabel('Power (W)', color='tab:blue')
        ax1.plot(wind_speeds, power_curve, color='tab:blue')
        ax1.tick_params(axis='y', labelcolor='tab:blue')
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.set_ylabel('Thrust (N)', color='tab:red')
        ax2.plot(wind_speeds, thrust_curve, color='tab:red')
        ax2.tick_params(axis='y', labelcolor='tab:red')

        plt.title('Power and Thrust Curves vs Wind Speed')
        plt.tight_layout()
        plt.show()

    def plot_pitch_rot_speed(self, wind_speeds):
        """
        Plot optimal pitch angle and rotational speed curves.

        Args:
            wind_speeds (np.ndarray): Array of wind speeds in m/s.
        """
        pitch_vals, rpm_vals = [], []

        for V0 in wind_speeds:
            pitch, rpm = self.get_optimal_operational_values(V0)
            pitch_vals.append(pitch)
            rpm_vals.append(rpm)

        plt.figure(figsize=(10, 4))
        plt.plot(wind_speeds, pitch_vals, label='Pitch Angle (deg)', color='green')
        plt.plot(wind_speeds, rpm_vals, label='Rotational Speed (rpm)',
                 color='purple')
        plt.xlabel('Wind Speed (m/s)')
        plt.ylabel('Pitch / RPM')
        plt.title('Optimal Pitch Angle and Rotational Speed vs Wind Speed')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()


