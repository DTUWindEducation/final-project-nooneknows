# Import the necessary libraries for numerical operations, data handling, and plotting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Define the BEMSolver class to handle Blade Element Momentum (BEM) calculations
class BEMSolver:
    def __init__(self, blade_data, airfoil_data, operational_data, B=3, rho=1.225):
        # Initialize blade, airfoil, and operational data
        self.blade_data = blade_data
        self.airfoil_data = airfoil_data
        self.operational_data = operational_data
        
        # Number of blades and air density
        self.B = B
        self.rho = rho
        
        # Maximum blade span and rotor swept area
        self.R = blade_data['BlSpn'].max()
        self.A = np.pi * self.R**2

    # Compute the tip speed ratio (TSR) given wind speed and rotational speed
    def compute_tsr(self, wind_speed, rot_speed_rpm):
        # Convert RPM to rad/s
        omega = (rot_speed_rpm * 2 * np.pi) / 60  
        # TSR formula
        return omega * self.R / wind_speed  

    # Compute lift and drag coefficients as a function of span position and angle of attack
    def compute_aero_coeff(self, r, alpha):
        # Find the closest index in blade data for the given span position
        idx = (np.abs(self.blade_data['BlSpn'] - r)).argmin()
        
        # Retrieve airfoil data for the corresponding airfoil ID
        airfoil_id = str(int(self.blade_data['BlAFID'][idx]) - 1).zfill(2)
        alpha_data = self.airfoil_data[airfoil_id]['alpha_deg']
        Cl_data = self.airfoil_data[airfoil_id]['Cl']
        Cd_data = self.airfoil_data[airfoil_id]['Cd']
        
        # Interpolate lift and drag coefficients based on the angle of attack
        Cl = np.interp(alpha, alpha_data, Cl_data)
        Cd = np.interp(alpha, alpha_data, Cd_data)
        return Cl, Cd

    # Compute axial (a) and tangential (a') induction factors
    def compute_induction(self, r, V0, theta_p, omega, c_r, beta_r, airfoil_id):
        # Initialize induction factors
        a, a_prime = 0, 0  
        # Blade solidity
        solidity = (c_r * self.B) / (2 * np.pi * r + 1e-6)  

        # Iteratively solve for induction factors
        for _ in range(100):
           # Flow angle 
            phi = np.arctan((1 - a) * V0 / ((1 + a_prime) * omega * r + 1e-6))  
            # Angle of attack
            alpha = np.degrees(phi) - (theta_p + np.rad2deg(beta_r))  
            # Aerodynamic coefficients
            Cl, Cd = self.compute_aero_coeff(r, alpha)  
            # Normal force coefficient
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)  
            # Tangential force coefficient
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)  
            
            # Update axial and tangential induction factors
            a_new = 1 / (4 * np.sin(phi)**2 / (solidity * Cn) + 1)
            a_prime_new = 1 / (4 * np.sin(phi) * np.cos(phi) / (solidity * Ct) - 1)
            
            # Check for convergence
            if np.abs(a - a_new) < 1e-4 and np.abs(a_prime - a_prime_new) < 1e-4:
                break
            a, a_prime = a_new, a_prime_new
        
        return a_new, a_prime_new

    # Solve BEM equations to compute thrust, torque, and power
    def solve_bem(self, V0, theta_p, rot_speed_rpm):
        # Convert RPM to rad/s
        omega = (rot_speed_rpm * 2 * np.pi) / 60  
        # Initialize thrust and torque
        thrust, torque = 0, 0  
        # Rated power in watts
        rated_power = 15678.725250 * 1e3  

        # Loop through each blade span section
        for i, r in enumerate(self.blade_data['BlSpn']):
            # Chord length
            c_r = self.blade_data['BlChord'][i]  
            # Twist angle in radians
            beta_r = np.radians(self.blade_data['BlTwist'][i]) 
            # Airfoil ID
            airfoil_id = str(int(self.blade_data['BlAFID'][i]) - 1).zfill(2)  
            
            # Compute induction factors
            a, a_prime = self.compute_induction(r, V0, theta_p, omega, c_r, beta_r, airfoil_id)
            # Span section length
            dr = self.R / len(self.blade_data)  
            
            # Compute differential thrust and torque
            dT = 4 * np.pi * r * self.rho * V0**2 * a * (1 - a) * dr
            dM = 4 * np.pi * r**3 * self.rho * V0 * omega * a_prime * (1 - a) * dr
            thrust += dT
            torque += dM
        
        # Compute power, thrust coefficient, and power coefficient
        power = torque * omega
        CT = thrust / (0.5 * self.rho * self.A * V0**2)
        CP = power / (0.5 * self.rho * self.A * V0**3)
        
        # Cap power at rated wind speed
        if V0 >= self.get_rated_wind_speed():
            power = rated_power
        
        return thrust, torque, power, CT, CP

    # Get optimal pitch angle and rotational speed for a given wind speed
    def get_optimal_operational_values(self, wind_speed):
         # Wind speed data
        wind_speeds = self.operational_data['wind_speed_ms'] 
        # Optimal pitch angles
        optimal_pitch_angles = self.operational_data['pitch_deg']  
        # Rotational speeds
        optimal_rot_speeds = self.operational_data['rot_speed_rpm']  
        # Maximum rotational speed
        max_rot_speed = np.max(optimal_rot_speeds)  

        # Find the rated wind speed and power
        indices = np.where(optimal_rot_speeds == max_rot_speed)
        index = indices[0][0]
        rated_wind_speed = wind_speeds[index]
        rated_power = self.operational_data['aero_power_kw'][index]

        # Interpolate pitch and rotational speed
        optimal_pitch = np.interp(wind_speed, wind_speeds, optimal_pitch_angles)
        optimal_rot_speed = np.interp(wind_speed, wind_speeds, optimal_rot_speeds)

        # Adjust values for wind speeds above rated wind speed
        if wind_speed >= rated_wind_speed:
            optimal_rot_speed = max_rot_speed
            optimal_pitch = np.interp(wind_speed, wind_speeds, optimal_pitch_angles)
        
        return optimal_pitch, optimal_rot_speed

    # Get the rated wind speed based on maximum rotational speed
    def get_rated_wind_speed(self):
        max_rot_speed = np.max(self.operational_data['rot_speed_rpm'])
        rated_index = np.where(self.operational_data['rot_speed_rpm'] == max_rot_speed)[0][0]
        return self.operational_data['wind_speed_ms'][rated_index]

    # Compute power and thrust curves for a range of wind speeds
    def compute_power_thrust_curve(self, wind_speeds):
        # Initialize power curve
        power_curve = []  
        # Initialize thrust curve
        thrust_curve = []  

        # Loop through wind speeds and compute power and thrust
        for V0 in wind_speeds:
            theta_p, omega_rpm = self.get_optimal_operational_values(V0)
            T, M, P, CT, CP = self.solve_bem(V0, theta_p, omega_rpm)
            power_curve.append(P)
            thrust_curve.append(T)
        
        return np.array(power_curve), np.array(thrust_curve)

    # Plot power and thrust curves as a function of wind speed
    def plot_power_thrust_curves(self, wind_speeds, power_curve, thrust_curve):
        fig, ax1 = plt.subplots()

        # Plot power curve on the primary y-axis
        ax1.set_xlabel('Wind Speed (m/s)')
        ax1.set_ylabel('Power (W)', color='tab:blue')
        ax1.plot(wind_speeds, power_curve, color='tab:blue')
        ax1.tick_params(axis='y', labelcolor='tab:blue')

        # Plot thrust curve on the secondary y-axis
        ax2 = ax1.twinx()
        ax2.set_ylabel('Thrust (N)', color='tab:red')
        ax2.plot(wind_speeds, thrust_curve, color='tab:red')
        ax2.tick_params(axis='y', labelcolor='tab:red')

        plt.title('Power and Thrust Curves vs Wind Speed')
        plt.tight_layout()
        plt.show()

    # Compute and plot optimal pitch angle and rotational speed curves
    def plot_pitch_rot_speed(self, wind_speeds):
        # Initialize pitch angle values
        pitch_vals = []  
        # Initialize rotational speed values
        rpm_vals = [] 

        # Loop through wind speeds and compute pitch and rotational speed
        for V0 in wind_speeds:
            pitch, rpm = self.get_optimal_operational_values(V0)
            pitch_vals.append(pitch)
            rpm_vals.append(rpm)

        # Plot the results
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


