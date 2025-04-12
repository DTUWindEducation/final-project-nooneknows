# Import classes
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
import numpy as np

import load_files, solve_bem , postprocessing

# import (
#     BladeGeometryLoader,
#     AirfoilDataLoader,
#     OperationalDataLoader,
#     BEMSolver,
#     plot_all_airfoils,
#     plot_power_thrust_curves,
# )

blade_loader = load_files.BladeGeometryLoader()
blade_data, _ = blade_loader.load()
airfoil_loader = load_files.AirfoilDataLoader()
airfoil_data = airfoil_loader.load()
operational_loader = load_files.OperationalDataLoader()
operational_data = operational_loader.load()
bem_solver = solve_bem.BEMSolver(blade_data, airfoil_data, operational_data)
thrust, torque, power, CT, CP = bem_solver.solve_bem(10, 0, 5)
# bem_solver.plot_power_thrust_curves()
postprocessing.plot_all_airfoils(airfoil_data)

wind_speed = 10  # m/s
pitch_angle = 0  # degrees
rot_speed_rpm = 5  # rpm

thrust, torque, power, CT, CP = bem_solver.solve_bem(wind_speed, pitch_angle, rot_speed_rpm)

print(f"Thrust: {thrust:.2f} N")
print(f"Torque: {torque:.2f} Nm")
print(f"Power: {power:.2f} W")
print(f"Thrust Coefficient (CT): {CT:.4f}")
print(f"Power Coefficient (CP): {CP:.4f}")
postprocessing.plot_power_thrust_curves(operational_data)


wind_speeds = np.linspace(3, 25, 100)
power_curve, thrust_curve = bem_solver.compute_power_thrust_curve(wind_speeds)
bem_solver.plot_power_thrust_curves(wind_speeds, power_curve, thrust_curve)

