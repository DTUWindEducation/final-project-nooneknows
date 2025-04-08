# Import classes
from WT import (
    BladeGeometryLoader,
    AirfoilDataLoader,
    OperationalDataLoader,
    BEMSolver,
    plot_all_airfoils,
    plot_power_thrust_curves,
)

blade_loader = BladeGeometryLoader()
blade_data, _ = blade_loader.load()
airfoil_loader = AirfoilDataLoader()
airfoil_data = airfoil_loader.load()
operational_loader = OperationalDataLoader()
operational_data = operational_loader.load()
bem_solver = BEMSolver(blade_data, airfoil_data, operational_data)
thrust, torque, power, CT, CP = bem_solver.solve_bem(10, 0, 5)
# bem_solver.plot_power_thrust_curves()
plot_all_airfoils(airfoil_data)

wind_speed = 10  # m/s
pitch_angle = 0  # degrees
rot_speed_rpm = 5  # rpm

thrust, torque, power, CT, CP = bem_solver.solve_bem(wind_speed, pitch_angle, rot_speed_rpm)

print(f"Thrust: {thrust:.2f} N")
print(f"Torque: {torque:.2f} Nm")
print(f"Power: {power:.2f} W")
print(f"Thrust Coefficient (CT): {CT:.4f}")
print(f"Power Coefficient (CP): {CP:.4f}")
plot_power_thrust_curves(operational_data)