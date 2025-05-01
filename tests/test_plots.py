"""
Unit tests for the postprocessing module.

These tests validate the functionality of the plotting functions for airfoil
shapes, multiple airfoils, and power/thrust curves.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
from wt import postprocessing as pp


@pytest.fixture(autouse=True)
def prevent_show(monkeypatch):
    """
    Prevent plt.show() from blocking tests by replacing it with a no-op.
    """
    monkeypatch.setattr(plt, "show", lambda: None)


def test_plot_airfoil_shape():
    """
    Test the plot_airfoil_shape function with a simple sinusoidal airfoil.
    """
    x = np.linspace(0, 1, 50)
    y = np.sin(2 * np.pi * x) * 0.1

    pp.plot_airfoil_shape(
        x, y, airfoil_id="TestAirfoil", color="green", linestyle="--"
    )
    # If no exception is raised, the test passes.


def test_plot_all_airfoils():
    """
    Test the plot_all_airfoils function with data for two airfoils.
    """
    airfoil_data = {
        "001": {
            "x_coords": np.linspace(0, 1, 50),
            "y_coords": np.sin(2 * np.pi * np.linspace(0, 1, 50)) * 0.1,
        },
        "002": {
            "x_coords": np.linspace(0, 1, 50),
            "y_coords": np.cos(2 * np.pi * np.linspace(0, 1, 50)) * 0.1,
        },
    }

    pp.plot_all_airfoils(airfoil_data)
    # If no exception is raised, the test passes.


def test_plot_power_thrust_curves():
    """
    Test the plot_power_thrust_curves function
    to ensure it runs with no errors
    with sample operational data.
    """
    operational_data = {
        "wind_speed_ms": np.array([3, 5, 8, 12, 15, 20]),
        "pitch_deg": np.array([2, 1.5, 1, 0.5, 0.2, 0]),
        "rot_speed_rpm": np.array([10, 20, 30, 40, 50, 60]),
        "aero_power_kw": np.array([0, 100, 500, 1000, 1500, 2000]),
        "aero_thrust_kw": np.array([1, 2, 3, 4, 5, 6]),
    }

    pp.plot_power_thrust_curves(operational_data)
    # If no exception is raised, the test passes.


def test_plot_spanwise_normal_tangential_loads():
    """
    Test the plot_spanwise_normal_tangential_loads function
     with sample data.
     """
    data = {
        "r": np.linspace(0, 50, 100),
        "Fn": np.linspace(1000, 2000, 100),
        "Ft": np.linspace(500, 1000, 100),
    }
    pp.plot_spanwise_normal_tangential_loads(data)
    plt.close()


def test_plot_pitch_rot_speed():
    """
    Test that the plot_pitch_rot_speed function
    runs without crashing and produces expected plot
    with sample data.
    """
    wind_speeds = np.linspace(4, 25, 50)
    pitch_vals = np.linspace(0, 20, 50)
    rpm_vals = np.linspace(10, 100, 50)

    pp.plot_pitch_rot_speed(wind_speeds, pitch_vals, rpm_vals)
    plt.close()


def test_plot_bem_power_thrust_curves():
    """
    Test that the plot_bem_power_thrust_curves function 
    execute without errors with sample data.
    """
    wind_speeds = np.linspace(4, 25, 50)
    power_curve = np.random.uniform(1e4, 1e6, 50)
    thrust_curve = np.random.uniform(1e3, 1e5, 50)

    pp.plot_bem_power_thrust_curves(wind_speeds, power_curve, thrust_curve)
    plt.close()
