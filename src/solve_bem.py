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
    
    def compute_aero_coeff(radius,alpha,alpha_data,Cl_data,Cd_data):
        Cl = np.interp(alpha, alpha_data, Cl_data)
        Cd = np.interp(alpha, alpha_data, Cd_data)
        return Cl, Cd
    
    def compute_induction(r, wind_speed, pitch_angle, omega, solidity, beta_r, alpha_data,Cl_data,Cd_data):
        a, a_prime = 0, 0
        for _ in range(100):
            phi = np.arctan((1-a) * wind_speed / ((1+a_prime) * omega * r+ 1e-6))
            alpha = np.degrees(phi) - (pitch_angle + np.rad2deg(beta_r))
            Cl, Cd = BEMSolver.compute_aero_coeff(r,alpha,alpha_data,Cl_data,Cd_data)
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
            a_new = 1 / (4 * np.sin(phi)**2 / (solidity * Cn) + 1)
            a_prime_new = 1 / (4 * np.sin(phi) * np.cos(phi) / (solidity * Ct) - 1)
            if np.abs(a - a_new) < 1e-4 and np.abs(a_prime - a_prime_new) < 1e-4:
                break
            a, a_prime = a_new, a_prime_new
        return a_new, a_prime_new

        

    def solve_bem(self, wind_speed, pitch_angle, rot_speed_rpm):
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        thrust, torque = 0, 0
        induction_factors = []  # store induction factors (a, a') as a 
        lift_coeff = []
        drag_coeff = []
        rated_power = 15678.725250*1e3 # To store rated power
        for i, r in enumerate(self.blade_data['BlSpn']):
            c_r = self.blade_data['BlChord'][i]
            beta_r = np.radians(self.blade_data['BlTwist'][i])
            airfoil_id = str(int(self.blade_data['BlAFID'][i]) - 1).zfill(2)
            alpha_data = self.airfoil_data[airfoil_id]['alpha_deg']
            Cl_data = self.airfoil_data[airfoil_id]['Cl']
            Cd_data = self.airfoil_data[airfoil_id]['Cd']
            solidity = (c_r * self.B) / (2 * np.pi * r+ 1e-6)
            a, a_prime = 0, 0
            
            a, a_prime = BEMSolver.compute_induction(r, wind_speed, pitch_angle, omega, solidity, beta_r, alpha_data,Cl_data,Cd_data)
            induction_factors.append((a, a_prime))  # Store induction factors for each span
            dT = 4 * np.pi * r * self.rho * wind_speed**2 * a * (1-a) * (self.R/len(self.blade_data))
            dM = 4 * np.pi * r**3 * self.rho * wind_speed * omega * a_prime * (1-a) * (self.R/len(self.blade_data))
            thrust += dT
            torque += dM
        power = torque * omega
        CT = thrust / (0.5 * self.rho * self.A * wind_speed**2)
        CP = power / (0.5 * self.rho * self.A * wind_speed**3)
        # If power exceeds rated power, cap it
        if rated_power is not None and wind_speed >= self.operational_data['wind_speed_ms'][np.argmax(self.operational_data['rot_speed_rpm'])]:
            power = rated_power

        return thrust, torque, power, CT, CP

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
    
     # 3. Function to compute power and thrust curve
    def compute_power_thrust_curve(self, wind_speeds):
        power_curve = []
        thrust_curve = []
        for wind_speed in wind_speeds:
            pitch_angle, rot_speed_rpm = self.get_optimal_operational_values(wind_speed)
            thrust, torque, power, CT, CP = self.solve_bem(wind_speed, pitch_angle, rot_speed_rpm)
            power_curve.append(power)
            thrust_curve.append(thrust)

        return np.array(power_curve), np.array(thrust_curve)

    # 4. Function to plot power and thrust curves
    def plot_power_thrust_curves(self, wind_speeds, power_curve, thrust_curve):
        fig, ax1 = plt.subplots()

        ax1.set_xlabel('Wind Speed (m/s)')
        ax1.set_ylabel('Power (W)', color='tab:blue')
        ax1.plot(wind_speeds, power_curve, color='tab:blue', label='Power Curve')
        ax1.tick_params(axis='y', labelcolor='tab:blue')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Thrust (N)', color='tab:red')
        ax2.plot(wind_speeds, thrust_curve, color='tab:red', label='Thrust Curve')
        ax2.tick_params(axis='y', labelcolor='tab:red')

        fig.tight_layout()
        plt.title('Power and Thrust Curves')
        plt.show()


