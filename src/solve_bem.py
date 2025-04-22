import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class BEMSolver:
    def __init__(self, blade_data, airfoil_data, operational_data, B=3, rho=1.225):
        self.blade_data = blade_data
        self.airfoil_data = airfoil_data
        self.operational_data = operational_data
        self.B = B
        self.rho = rho
        self.R = blade_data['BlSpn'].max()
        self.A = np.pi * self.R**2

    def compute_tsr(self, wind_speed, rot_speed_rpm):
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        return omega * self.R / wind_speed
    
    # 3. Lift and drag coefficient as function of r and alpha
    def compute_aero_coeff(self, r, alpha):
        # Find the closest index in blade_data for the given span position r
        idx = (np.abs(self.blade_data['BlSpn'] - r)).argmin()
        
        # Lookup corresponding airfoil ID and get data
        airfoil_id = str(int(self.blade_data['BlAFID'][idx]) - 1).zfill(2)
        alpha_data = self.airfoil_data[airfoil_id]['alpha_deg']
        Cl_data = self.airfoil_data[airfoil_id]['Cl']
        Cd_data = self.airfoil_data[airfoil_id]['Cd']
        
        # Interpolate aerodynamic coefficients
        Cl = np.interp(alpha, alpha_data, Cl_data)
        Cd = np.interp(alpha, alpha_data, Cd_data)
        return Cl, Cd

    
    # 4. Compute a and a' as function of r, V0, θp, ω
    def compute_induction(self, r, V0, theta_p, omega, c_r, beta_r, airfoil_id):
        a, a_prime = 0, 0
        solidity = (c_r * self.B) / (2 * np.pi * r + 1e-6)
        for _ in range(100):
            phi = np.arctan((1 - a) * V0 / ((1 + a_prime) * omega * r + 1e-6))
            alpha = np.degrees(phi) - (theta_p + np.rad2deg(beta_r))
            Cl, Cd = self.compute_aero_coeff(r, alpha)
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
            a_new = 1 / (4 * np.sin(phi)**2 / (solidity * Cn) + 1)
            a_prime_new = 1 / (4 * np.sin(phi) * np.cos(phi) / (solidity * Ct) - 1)
            if np.abs(a - a_new) < 1e-4 and np.abs(a_prime - a_prime_new) < 1e-4:
                break
            a, a_prime = a_new, a_prime_new
        return a_new, a_prime_new

        

     # 5. Compute T, M, P as function of V0, θp, ω
    def solve_bem(self, V0, theta_p, rot_speed_rpm):
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        thrust, torque = 0, 0
        rated_power = 15678.725250 * 1e3
        for i, r in enumerate(self.blade_data['BlSpn']):
            c_r = self.blade_data['BlChord'][i]
            beta_r = np.radians(self.blade_data['BlTwist'][i])
            airfoil_id = str(int(self.blade_data['BlAFID'][i]) - 1).zfill(2)
            a, a_prime = self.compute_induction(r, V0, theta_p, omega, c_r, beta_r, airfoil_id)
            dr = self.R / len(self.blade_data)
            dT = 4 * np.pi * r * self.rho * V0**2 * a * (1 - a) * dr
            dM = 4 * np.pi * r**3 * self.rho * V0 * omega * a_prime * (1 - a) * dr
            thrust += dT
            torque += dM
        power = torque * omega
        CT = thrust / (0.5 * self.rho * self.A * V0**2)
        CP = power / (0.5 * self.rho * self.A * V0**3)
        if V0 >= self.get_rated_wind_speed():
            power = rated_power
        return thrust, torque, power, CT, CP


    # 6. Compute θp and ω as function of V0 using strategy
    def get_optimal_operational_values(self, wind_speed):
        # Assuming the operational strategy contains wind speed, optimal pitch angle, and rotational speed
        wind_speeds = self.operational_data['wind_speed_ms']
        optimal_pitch_angles = self.operational_data['pitch_deg']
        optimal_rot_speeds = self.operational_data['rot_speed_rpm']
        max_rot_speed = np.max(optimal_rot_speeds)

        # Find the indices where the value is max rotational speed
        indices = np.where(optimal_rot_speeds == max_rot_speed)

        # Get the first occurrence index (if there are multiple)
        index = indices[0][0]

        # Use np.interp for linear interpolation
        optimal_pitch = np.interp(wind_speed, wind_speeds, optimal_pitch_angles)
        optimal_rot_speed = np.interp(wind_speed, wind_speeds, optimal_rot_speeds)

        # If the wind speed is above or equal to the rated wind speed, 
        # ensure power remains constant and cap the rotational speed and power.
        rated_wind_speed = wind_speeds[index]  # Rated wind speed at which max rotation speed occurs
        rated_power = self.operational_data['aero_power_kw'][index]  # Rated power at that wind speed

        if wind_speed >= rated_wind_speed:
            # After rated wind speed, power remains constant, so we do not change rotational speed
            optimal_rot_speed = max_rot_speed
            optimal_pitch = np.interp(wind_speed, wind_speeds, optimal_pitch_angles)  # Still adjust pitch
        else:
            # Before rated wind speed, calculate optimal pitch and rotational speed as normal
            optimal_pitch = np.interp(wind_speed, wind_speeds, optimal_pitch_angles)
            optimal_rot_speed = np.interp(wind_speed, wind_speeds, optimal_rot_speeds)

        return optimal_pitch, optimal_rot_speed
    
    def get_rated_wind_speed(self):
        max_rot_speed = np.max(self.operational_data['rot_speed_rpm'])
        rated_index = np.where(self.operational_data['rot_speed_rpm'] == max_rot_speed)[0][0]
        return self.operational_data['wind_speed_ms'][rated_index]

    
    # 7. Compute power and thrust curve P(V0), T(V0)
    def compute_power_thrust_curve(self, wind_speeds):
        power_curve = []
        thrust_curve = []
        for V0 in wind_speeds:
            theta_p, omega_rpm = self.get_optimal_operational_values(V0)
            T, M, P, CT, CP = self.solve_bem(V0, theta_p, omega_rpm)
            power_curve.append(P)
            thrust_curve.append(T)
        return np.array(power_curve), np.array(thrust_curve)

    def plot_power_thrust_curves(self, wind_speeds, power_curve, thrust_curve):
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Wind Speed (m/s)')
        ax1.set_ylabel('Power (W)', color='tab:blue')
        ax1.plot(wind_speeds, power_curve, color='tab:blue')
        ax1.tick_params(axis='y', labelcolor='tab:blue')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Thrust (N)', color='tab:red')
        ax2.plot(wind_speeds, thrust_curve, color='tab:red')
        ax2.tick_params(axis='y', labelcolor='tab:red')

        plt.title('Power and Thrust Curves vs Wind Speed')
        plt.tight_layout()
        plt.show()

    # === Extra Utility Function 1 ===
    def compute_tip_speed_ratio(self, wind_speed, rot_speed_rpm):
        return self.compute_tsr(wind_speed, rot_speed_rpm)

    def compute_tsr(self, wind_speed, rot_speed_rpm):
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        return omega * self.R / wind_speed

    # === Extra Utility Function 2 ===
    def plot_pitch_rot_speed(self, wind_speeds):
        pitch_vals = []
        rpm_vals = []
        for V0 in wind_speeds:
            pitch, rpm = self.get_optimal_operational_values(V0)
            pitch_vals.append(pitch)
            rpm_vals.append(rpm)

        plt.figure(figsize=(10, 4))
        plt.plot(wind_speeds, pitch_vals, label='Pitch Angle (deg)', color='green')
        plt.plot(wind_speeds, rpm_vals, label='Rotational Speed (rpm)', color='purple')
        plt.xlabel('Wind Speed (m/s)')
        plt.ylabel('Pitch / RPM')
        plt.title('Optimal Pitch Angle and Rotational Speed vs Wind Speed')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()


