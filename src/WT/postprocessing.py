import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import numpy as np
import pandas as pd


def plot_airfoil_shape(x_coords, y_coords, airfoil_id="Airfoil", color="blue", linestyle="-"):
    """Plot the airfoil shape."""
    plt.plot(x_coords, y_coords, label=airfoil_id, color=color, linestyle=linestyle, linewidth=2)

def plot_all_airfoils(airfoil_data):
    """Plot the shape of all airfoils on the same figure."""
    plt.figure(figsize=(12, 8))  # Create a single figure with larger dimensions

    # Generate a colormap to apply unique colors
    cmap = cm.get_cmap('tab20')  # Use 'tab20' colormap (you can change to other colormaps)
    
    # Plot each airfoil shape
    for i, (af_key, airfoil) in enumerate(airfoil_data.items()):
        color = cmap(i % 20)  # Cycle through colors
        linestyle = "--" if i % 2 == 0 else "-"  # Alternate line styles for variety
        plot_airfoil_shape(airfoil['x_coords'], airfoil['y_coords'], airfoil_id=f"Airfoil {af_key}", color=color, linestyle=linestyle)
    
    # Customize the plot
    plt.title("Airfoil Shapes", fontsize=18, fontweight='bold')
    plt.xlabel("x-coordinate", fontsize=14)
    plt.ylabel("y-coordinate", fontsize=14)
    plt.gca().set_aspect('equal', adjustable='box')  # Make the aspect ratio equal
    plt.grid(True, linestyle="--", alpha=0.7)  # Subtle grid lines
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # plt.legend(loc="upper left", fontsize=12, bbox_to_anchor=(1.05, 1))  # Place legend outside the plot
    
    plt.tight_layout()  # Ensure the plot fits within the figure area
    plt.show()

def plot_power_thrust_curves(operational_data):
        wind_speeds = np.arange(0.5, 25.5, 0.1)  # Wind speeds from 0.5 to 25 m/s with step 0.5
        
        # Extract operational data from self.operational_data
        wind_speed_ms = operational_data['wind_speed_ms']
        pitch_deg = operational_data['pitch_deg']
        rot_speed_rpm = operational_data['rot_speed_rpm']
        aero_power_kw = operational_data['aero_power_kw']
        aero_thrust_kn = operational_data['aero_thrust_kw']  # The provided data uses 'kw' but should be 'kn' for thrust
        
        # Interpolate pitch, rotational speed, power, and thrust
        pitch_interp = np.interp(wind_speeds, wind_speed_ms, pitch_deg)
        omega_interp = np.interp(wind_speeds, wind_speed_ms, rot_speed_rpm)
        power_interp = np.interp(wind_speeds, wind_speed_ms, aero_power_kw)
        thrust_interp = np.interp(wind_speeds, wind_speed_ms, aero_thrust_kn)

        # Create side-by-side subplots
        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        # Plot Power Curve
        axs[0].plot(wind_speeds, power_interp, label='Power (kW)', color='b')
        axs[0].set_xlabel('Wind Speed (m/s)')
        axs[0].set_ylabel('Power (kW)')
        axs[0].set_title('Power Curve')
        axs[0].grid(True)
        axs[0].legend()

        # Plot Thrust Curve
        axs[1].plot(wind_speeds, thrust_interp, label='Thrust (kN)', color='r')
        axs[1].set_xlabel('Wind Speed (m/s)')
        axs[1].set_ylabel('Thrust (kN)')
        axs[1].set_title('Thrust Curve')
        axs[1].grid(True)
        axs[1].legend()

        # Adjust layout
        plt.tight_layout()
        plt.show()