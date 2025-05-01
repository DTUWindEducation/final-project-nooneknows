"""
Module for Blade Element Momentum (BEM) calculations.

This module defines the BEMSolver class, which provides methods for computing
aerodynamic coefficients, induction factors, and power/thrust curves for wind
turbines.
"""

import numpy as np
# pylint: disable=C0103


class BEMSolver:
    """
    A class to handle Blade Element Momentum (BEM) calculations.
    """

    def __init__(self, blade_data, airfoil_data,
                 operational_data, B=3, rho=1.225):
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
        cl_data = self.airfoil_data[airfoil_id]['Cl']
        cd_data = self.airfoil_data[airfoil_id]['Cd']
        cl = np.interp(alpha, alpha_data, cl_data)
        cd = np.interp(alpha, alpha_data, cd_data)
        return cl, cd

    def compute_induction(self, r, wind_speed, pitch_angle,
                          omega, chord, twist, airfoil_id):
        """
        Compute axial (a) and tangential (a') induction factors.

        Args:
            r (float): Span position in meters.
            wind_speed (float): Wind speed in m/s.
            pitch_angle (float): Pitch angle in degrees.
            omega (float): Rotational speed in rad/s.
            chord (float): Chord length in meters.
            twist (float): Twist angle in radians.
            airfoil_id (str): Airfoil ID.

        Returns:
            tuple: Axial induction factor (a)
              and tangential induction factor (a').
        """
        a, a_prime = 0, 0
        solidity = (chord * self.B) / (2 * np.pi * r + 1e-6)

        for _ in range(100):
            phi = np.arctan((1 - a) * wind_speed /
                            ((1 + a_prime) * omega * r + 1e-6))
            alpha = np.degrees(phi) - (pitch_angle + np.rad2deg(twist))
            cl, cd = self.compute_aero_coeff(r, alpha)
            cn = cl * np.cos(phi) + cd * np.sin(phi)
            ct = cl * np.sin(phi) - cd * np.cos(phi)
            a_new = 1 / (4 * np.sin(phi)**2 / (solidity * cn) + 1)
            a_prime_new = 1 / (4 * np.sin(phi) * np.cos(phi) /
                               (solidity * ct) - 1)
            if (np.abs(a - a_new) < 1e-4 and
                    np.abs(a_prime - a_prime_new) < 1e-4):
                break
            a, a_prime = a_new, a_prime_new

        return a_new, a_prime_new

    def solve_bem(self, wind_speed, pitch_angle, rot_speed_rpm,
                  target_power_kw=None):
        """
        Solve BEM equations to compute thrust, torque, and power.

        Args:
            wind_speed (float): Wind speed in m/s.
            pitch_angle (float): Pitch angle in degrees.
            rot_speed_rpm (float): Rotational speed in RPM.
            target_power_kw (float, optional): Target power in kW.

        Returns:
            tuple: Thrust (N), torque (Nm), power (W), thrust coefficient (CT),
                   and power coefficient (CP).
        """
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        thrust, torque = 0, 0
        dr = self.R / len(self.blade_data['BlSpn'])
        dT_vals, dM_vals = [], []

        for i, r in enumerate(self.blade_data['BlSpn']):
            chord = self.blade_data['BlChord'][i]
            twist = np.radians(self.blade_data['BlTwist'][i])
            airfoil_id = str(int(self.blade_data['BlAFID'][i]) - 1).zfill(2)
            a, a_prime = self.compute_induction(r, wind_speed, pitch_angle,
                                                omega, chord, twist, airfoil_id)
            dT = 4 * np.pi * r * self.rho * wind_speed**2 * a * (1 - a) * dr
            dM = (4 * np.pi * r**3 * self.rho * wind_speed * omega * a_prime *
                  (1 - a) * dr)
            dT_vals.append(dT)
            dM_vals.append(dM)

        thrust = np.trapz(dT_vals, dx=dr)
        torque = np.trapz(dM_vals, dx=dr)
        power = torque * omega

        if target_power_kw is not None:
            target_power_kw = target_power_kw * 1e3
            power = min(power, target_power_kw)

        ct = thrust / (0.5 * self.rho * self.A * wind_speed**2)
        cp = power / (0.5 * self.rho * self.A * wind_speed**3)

        return thrust, torque, power, ct, cp

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
        interp_pitch = np.interp(
            wind_speed, wind_speeds, optimal_pitch_angles)
        interp_rot_speed = np.interp(
            wind_speed, wind_speeds, optimal_rot_speeds)

        # Identify rated power and rated wind speed
        rated_power = np.max(power_curve_kw)
        rated_indices = np.where(power_curve_kw == rated_power)[0]
        rated_index = rated_indices[0]
        rated_wind_speed = wind_speeds[rated_index]
        rated_rot_speed = optimal_rot_speeds[rated_index]

        # If wind speed >= rated, cap at rated power behavior
        if wind_speed >= rated_wind_speed:
            # Pitch still varies to regulate power
            optimal_pitch = interp_pitch
            # Fix RPM at rated
            optimal_rot_speed = rated_rot_speed
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
            interp_power_kw = np.interp(V0,
                                        self.operational_data['wind_speed_ms'], power_curve_kw)
            T, _ , P, _, _ = self.solve_bem(V0, theta_p,
                                              omega_rpm, target_power_kw=interp_power_kw)
            power_curve.append(P)
            thrust_curve.append(T)
            omega_curve.append(omega_rpm)
            pitch_curve.append(theta_p)
        return np.array(power_curve), np.array(thrust_curve), np.array(omega_curve), np.array(pitch_curve)

    def compute_spanwise_normal_tangential_loads(self, V0, theta_p, rot_speed_rpm):
        """
        N.1 EXTRA FUNCTION:
        Compute spanwise normal and tangential force distributions (per unit length).

        Args:
            V0 (float): Wind speed in m/s.
            theta_p (float): Pitch angle in degrees.
            rot_speed_rpm (float): Rotational speed in RPM.

        Returns:
            dict: Dictionary with span positions and normal/tangential load components.
        """
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        r_vals = self.blade_data['BlSpn']
        normal_loads, tangential_loads = [], []

        for i, r in enumerate(r_vals):
            c_r = self.blade_data['BlChord'][i]
            beta_r = np.radians(self.blade_data['BlTwist'][i])
            airfoil_id = str(int(self.blade_data['BlAFID'][i]) - 1).zfill(2)

            a, a_prime = self.compute_induction(r, V0, theta_p, omega, c_r, beta_r, airfoil_id)

            phi = np.arctan((1 - a) * V0 / ((1 + a_prime) * omega * r + 1e-6))
            alpha = np.degrees(phi) - (theta_p + np.rad2deg(beta_r))
            Cl, Cd = self.compute_aero_coeff(r, alpha)
            W = np.sqrt(((1 - a) * V0)**2 + ((1 + a_prime) * omega * r)**2)

            L = 0.5 * self.rho * W**2 * c_r * Cl
            D = 0.5 * self.rho * W**2 * c_r * Cd

            Fn = L * np.cos(phi) + D * np.sin(phi)  # Normal force (N/m)
            Ft = L * np.sin(phi) - D * np.cos(phi)  # Tangential force (N/m)

            normal_loads.append(Fn)
            tangential_loads.append(Ft)

        return {
            "r": r_vals,
            "Fn": np.array(normal_loads),
            "Ft": np.array(tangential_loads)
        }
