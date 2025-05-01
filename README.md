[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/zjSXGKeR)
# Our Great Package

Team: NoOneKnows

## Overview

This project uses the BEM theory to model the wind turbine. It takes as inputs:
 - blade geometry data - ./inputs/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat
 - airfoil data - ./inputs/IEA-15-240-RWT/Airfoils
 - operational data - IEA-15-240-RWT_AeroDyn15_Polar_{af_key}.dat

The purpose of the code is to forecast the aerodynamic behaviour of a wind turbine blade, and the final output is to plot the power curve and thrust curve.


The code performs the following key tasks:
- Loads wind turbine data.
- Computes lift coefficients and drag coefficients.
- Computes the axial and tangential induction factors.
- Computes the thrust torque and power.
- Computes optimal operational strategy.
- Computes power and thrust.
- Compute spanwise loads (tangential and normal)
- Generates plots.

## Quick-start guide

### 1. Clone the repository
```sh
git clone <repository_url>
cd <repository_directory>
```

### 2. Install the package
Create a new environment through the file environment.yml
```sh
conda env create -f environment.yml
```
Activate it
```sh
conda activate wt
```

Open the terminal and install the package:
```sh
pip install .
```

### 3. Install packages
Ensure you have Python installed. Then, install the required packages (numpy, scipy, matplotlib ecc..):

```sh
pip install numpy pandas pyplot
```

### 4. Run the main script
Execute the script to analyze the wind data:
```sh
python main.py
```
## 5. How the code works

The project consists of the following core components:

### **Main Script**
- **`main.py`**: The primary script that manages the entire process using the helper functions contained in **`__init__.py`** with the following steps:
  - Loads wind data from the `./inputs/IEA-15-240-RWT/` directory.
  - Computes the power and thrust values.
  - Generates plots to visualize results.

### **Helper Functions**
- **`__init__.py`**: Contains the following features:

  ### Data Loading:
  - BladeGeometryLoader: loads the blade geometry data.
  - AirfoilDataLoader: loads the airfoil data. 
  - OperationalDataLoader: loads the operational data. 

  ### Simulation:
  - BEMSolver: this class computes the rotational speed as well as the thrust, torque, power, lift and drag coefficients.

  ### Results Processing & Visualization:
- `plot_airfoil_shape()` plots the airfoil shape
- `plot_all_airfoils()` plots the airfoil shapes for all airfoils
- `plot_power_thrust_curves()` plots the thrust and power curves from operational data
- `plot_pitch_rot_speed()` plots rotational speed and pitch settings
- `plot_spanwise_normal_tangential_loads()` plots normal and tangential loads versus blade span
- `plot_bem_power_thrust_curves()` plot power and thrust coming from the bem code

### **Project Directories & files**
- **`./examples/`**: 
    - **main** 
- **`./inputs/`**
- **`./outputs/`** 
- **`./src`**
    - **`./__init__/`**
    - **`./load_files/`**
    - **`./postprocessing/`**
    - **`./solve_bem/`**
- **`project_structure.drawio`**: diagram of the code structure.


## Architecture

The package takes the blade geometry data, the airfoil data and the operational data through three different classes and uses them in the SolveBEM class. These classes are contained in main, where also the plots for the power curve and thrust curve are created.  

![alt text](inputs/Diagram.png)

The diagram source file can be found [here](./project_structure.drawio).

## Peer review

For this project, everybody tried to understand and begin to approach the code alone. This to have a deep understanding of the problem and a structured idea on how to solve it. 
After this, we discussed how to complete the project and together we wrote both the files contained in the **`__init__.py`** and **`main.py`**. This way, everybody had a clear view of the code and were able to understand the solution.
This, README file as well as the other that are not afore-mentioned were also produced together as a group.
