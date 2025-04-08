import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

    def solve_bem(self, wind_speed, pitch_angle, rot_speed_rpm):
        omega = (rot_speed_rpm * 2 * np.pi) / 60
        thrust, torque = 0, 0
        for i, r in enumerate(self.blade_data['BlSpn']):
            c_r = self.blade_data['BlChord'][i]
            beta_r = np.radians(self.blade_data['BlTwist'][i])
            airfoil_id = str(int(self.blade_data['BlAFID'][i]) - 1).zfill(2)
            alpha_data = self.airfoil_data[airfoil_id]['alpha_deg']
            Cl_data = self.airfoil_data[airfoil_id]['Cl']
            Cd_data = self.airfoil_data[airfoil_id]['Cd']
            solidity = (c_r * self.B) / (2 * np.pi * r+ 1e-6)
            a, a_prime = 0, 0
            for _ in range(100):
                phi = np.arctan((1-a) * wind_speed / ((1+a_prime) * omega * r+ 1e-6))
                alpha = np.degrees(phi) - (pitch_angle + np.rad2deg(beta_r))
                Cl = np.interp(alpha, alpha_data, Cl_data)
                Cd = np.interp(alpha, alpha_data, Cd_data)
                Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
                Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
                a_new = 1 / (4 * np.sin(phi)**2 / (solidity * Cn) + 1)
                a_prime_new = 1 / (4 * np.sin(phi) * np.cos(phi) / (solidity * Ct) - 1)
                if np.abs(a - a_new) < 1e-4 and np.abs(a_prime - a_prime_new) < 1e-4:
                    break
                a, a_prime = a_new, a_prime_new
            dT = 4 * np.pi * r * self.rho * wind_speed**2 * a * (1-a) * (self.R/len(self.blade_data))
            dM = 4 * np.pi * r**3 * self.rho * wind_speed * omega * a_prime * (1-a) * (self.R/len(self.blade_data))
            thrust += dT
            torque += dM
        power = torque * omega
        CT = thrust / (0.5 * self.rho * self.A * wind_speed**2)
        CP = power / (0.5 * self.rho * self.A * wind_speed**3)
        return thrust, torque, power, CT, CP

    