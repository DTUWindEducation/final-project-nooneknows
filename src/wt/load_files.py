"""
Module for loading and processing blade geometry, airfoil, and operational data.
"""

from pathlib import Path
import numpy as np
import pandas as pd


class BladeGeometryLoader:
    """
    A class to load and process blade geometry data from a file.
    """

    def __init__(self, filepath):
        """
        Initialize the BladeGeometryLoader with the file path.

        Args:
            filepath (str): Path to the blade geometry file.
        """
        self.filepath = filepath

    def load(self):
        """
        Load and process the blade geometry data.

        Returns:
            tuple: A tuple containing:
                - opt_data (pd.DataFrame): DataFrame with the blade geometry data.
                - units_opt (dict): Dictionary mapping column names to their units.
        """
        with open(self.filepath, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            index = next(i for i, line in enumerate(lines) if "BlSpn" in line)
            header_line = lines[index].strip()
            names = header_line.split()
            units = lines[index + 1].strip().split()

        start_index = next(i for i, line in enumerate(lines) if "(m)" in line)
        data = np.loadtxt(self.filepath, skiprows=start_index + 1)
        opt_data = pd.DataFrame(data, columns=names)
        units_opt = {names[i]: units[i] for i in range(len(names))}

        return opt_data, units_opt


class AirfoilDataLoader:
    """
    A class to load and process airfoil data from a directory.
    """

    def __init__(self, directory):
        """
        Initialize the AirfoilDataLoader with the directory path.

        Args:
            directory (str): Path to the directory containing airfoil data files.
        """
        self.directory = Path(directory)

    def load(self):
        """
        Load and process airfoil data.

        Returns:
            dict: A dictionary where keys are airfoil IDs and values are
                  dictionaries containing airfoil data (alpha, Cl, Cd,
                  x_coords, y_coords).
        """
        airfoil_data = {}
        for af_id in range(50):
            af_key = f"{af_id:02d}"
            polar_file = (
                self.directory / f"IEA-15-240-RWT_AeroDyn15_Polar_{af_key}.dat"
            )

            with open(polar_file, 'r', encoding='utf-8') as file:
                lines = file.readlines()
            start_index = next(
                i for i, line in enumerate(lines) if "!    (deg)" in line.lower()
            ) + 2
            polar_data = np.loadtxt(polar_file, skiprows=start_index)

            alpha, Cl, Cd = polar_data[:, 0], polar_data[:, 1], polar_data[:, 2]
            coords_file = (
                self.directory / f"IEA-15-240-RWT_AF{af_key}_Coords.txt"
            )
            x_coords, y_coords = self.load_airfoil_coordinates(coords_file)

            airfoil_data[af_key] = {
                'alpha_deg': alpha,
                'Cl': Cl,
                'Cd': Cd,
                'x_coords': x_coords,
                'y_coords': y_coords,
            }
        return airfoil_data

    @staticmethod
    def load_airfoil_coordinates(file_path):
        """
        Load airfoil coordinates from a file.

        Args:
            file_path (str): Path to the airfoil coordinates file.

        Returns:
            tuple: A tuple containing:
                - x_coords (list): List of x-coordinates.
                - y_coords (list): List of y-coordinates.
        """
        x_coords, y_coords = [], []
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            start_index = next(
                i for i, line in enumerate(lines) if "!  x/c " in line
            ) + 2
            for line in lines[start_index:]:
                if line.strip():
                    parts = line.split()
                    if len(parts) == 2:
                        x_coords.append(float(parts[0]))
                        y_coords.append(float(parts[1]))
        return x_coords, y_coords


class OperationalDataLoader:
    """
    A class to load and process operational data from a file.
    """

    def __init__(self, filepath):
        """
        Initialize the OperationalDataLoader with the file path.

        Args:
            filepath (str): Path to the operational data file.
        """
        self.filepath = filepath

    def load(self):
        """
        Load and process operational data.

        Returns:
            dict: A dictionary containing operational data:
                - wind_speed_ms (np.ndarray): Wind speed in m/s.
                - pitch_deg (np.ndarray): Pitch angle in degrees.
                - rot_speed_rpm (np.ndarray): Rotor speed in RPM.
                - aero_power_kw (np.ndarray): Aerodynamic power in kW.
                - aero_thrust_kw (np.ndarray): Aerodynamic thrust in kW.
        """
        with open(self.filepath, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        start_index = next(
            i for i, line in enumerate(lines)
            if "wind speed [m/s]" in line.lower()
        ) + 1
        operational_data = np.loadtxt(self.filepath, skiprows=start_index)

        return {
            'wind_speed_ms': operational_data[:, 0],
            'pitch_deg': operational_data[:, 1],
            'rot_speed_rpm': operational_data[:, 2],
            'aero_power_kw': operational_data[:, 3],
            'aero_thrust_kw': operational_data[:, 4],
        }
