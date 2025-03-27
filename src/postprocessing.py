import matplotlib.pyplot as plt
import matplotlib.cm as cm

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