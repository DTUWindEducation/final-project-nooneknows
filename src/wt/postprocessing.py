"""
Module for postprocessing and visualization of airfoil shapes and operational
data, including power and thrust curves.
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_airfoil_shape(
    x_coords, y_coords, airfoil_id="Airfoil", color="blue", linestyle="-"
):
    """
    Plot the shape of a single airfoil.

    Args:
        x_coords (list): List of x-coordinates.
        y_coords (list): List of y-coordinates.
        airfoil_id (str): Identifier for the airfoil (default: "Airfoil").
        color (str): Line color for the plot (default: "blue").
        linestyle (str): Line style for the plot (default: "-").
    """
    plt.plot(
        x_coords,
        y_coords,
        label=airfoil_id,
        color=color,
        linestyle=linestyle,
        linewidth=2,
    )


def plot_all_airfoils(airfoil_data):
    """
    Plot the shapes of all airfoils on a single figure.

    Args:
        airfoil_data (dict): Dictionary containing airfoil data. Keys are
            airfoil IDs, and values are dictionaries with 'x_coords' and
            'y_coords'.
    """
    plt.figure(figsize=(12, 8))
    cmap = plt.colormaps.get_cmap("tab20")

    for i, (af_key, airfoil) in enumerate(airfoil_data.items()):
        color = cmap(i % 20)
        linestyle = "--" if i % 2 == 0 else "-"
        plot_airfoil_shape(
            airfoil["x_coords"],
            airfoil["y_coords"],
            airfoil_id=f"Airfoil {af_key}",
            color=color,
            linestyle=linestyle,
        )

    plt.title("Airfoil Shapes", fontsize=18, fontweight="bold")
    plt.xlabel("x-coordinate", fontsize=14)
    plt.ylabel("y-coordinate", fontsize=14)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.show()


def plot_power_thrust_curves(operational_data):
    """
    Plot power and thrust curves based on operational data.

    Args:
        operational_data (dict): Dictionary containing operational data with
            keys:
            - 'wind_speed_ms': Wind speed in m/s.
            - 'pitch_deg': Pitch angle in degrees.
            - 'rot_speed_rpm': Rotor speed in RPM.
            - 'aero_power_kw': Aerodynamic power in kW.
            - 'aero_thrust_kw': Aerodynamic thrust in kW.
    """
    wind_speeds = np.arange(0.5, 25.5, 0.1)

    wind_speed_ms = operational_data["wind_speed_ms"]
    pitch_deg = operational_data["pitch_deg"]
    rot_speed_rpm = operational_data["rot_speed_rpm"]
    aero_power_kw = operational_data["aero_power_kw"]
    aero_thrust_kn = operational_data["aero_thrust_kw"]

    pitch_interp = np.interp(wind_speeds, wind_speed_ms, pitch_deg)
    omega_interp = np.interp(wind_speeds, wind_speed_ms, rot_speed_rpm)
    power_interp = np.interp(wind_speeds, wind_speed_ms, aero_power_kw)
    thrust_interp = np.interp(wind_speeds, wind_speed_ms, aero_thrust_kn)

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

    axs[0].plot(wind_speeds, power_interp, label="Power (kW)", color="b")
    axs[0].set_xlabel("Wind Speed (m/s)")
    axs[0].set_ylabel("Power (kW)")
    axs[0].set_title("Power Curve")
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(wind_speeds, thrust_interp, label="Thrust (kN)", color="r")
    axs[1].set_xlabel("Wind Speed (m/s)")
    axs[1].set_ylabel("Thrust (kN)")
    axs[1].set_title("Thrust Curve")
    axs[1].grid(True)
    axs[1].legend()

    plt.tight_layout()
    plt.show()