"""
This script was written for the purpose of assigments for the Combustion course on
Technical Univerity Munich.

Assigment: Exercise 6 - Ignition Delay Calculation
Student: Herman Andersson

"""

import pandas as pd
import numpy as np
import time
import os

# needed to save the figures without displaying them
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

import matplotlib.pyplot as plt
import cantera as ct

# -------------------- Experimental Data import Section --------------------------------
# Import experimental flame speed data
script_dir = os.path.dirname(__file__)
save_path = script_dir  
# save_path = "/path/to/your/desired/folder" 

IDT_05 = pd.read_csv(os.path.join(script_dir, 'IDT_0.5.csv'))
IDT_1 = pd.read_csv(os.path.join(script_dir, 'IDT_1.0.csv'))
IDT_2 = pd.read_csv(os.path.join(script_dir, 'IDT_2.0.csv'))
exp_datasets = [IDT_05, IDT_1, IDT_2]
# -------------------- Cantera Calculation Section --------------------------------
# Simulation Parameters
phi = [0.5,1.0,2.0]  # equivalence ratio range
T = np.linspace(800, 2000, 10)  # temperature in K
P = [1, 3, 5, 10, 25, 44]  # pressure in atm
P_pa = [p * ct.one_atm for p in P]  # Convert to Pa using Cantera's conversion
xaxis = 1000/T  # x-axis for the plot, 1000/T in K^-1


estimated_IDT = 1.0  # [s] Estimated ignition delay time in seconds
fuel = 'CH4:1'  # Fuel for the combustion reaction
gas = ct.Solution("gri30.yaml")  # Solution object used to compute mixture properties
ref_species = 'CH'  # Reference species for IDT determination

# Empty arrays to store results
store_arrays = []
for _ in range(3):
    store_arrays.append(np.zeros((len(P), len(T))))


def calculate_Ignition_Delay(pressure, temperature, gas, phi, IDT):
    for i in range(len(pressure)):
        for j in range(len(temperature)):
            print(f"Calculating IDT for P={pressure[i]/101325:.1f} atm, T={temperature[j]} K, phi={phi}")
            
            # Create oxidiser string for phi
            O2_str = str(2 / phi)
            N2_str = str(7.52 / phi)
            oxidiser = f'O2:{O2_str}, N2:{N2_str}'
            
            gas.set_equivalence_ratio(phi, fuel, oxidiser)
            gas.TP = temperature[j], pressure[i]
            
            reactor = ct.Reactor(contents=gas)
            reactor_network = ct.ReactorNet([reactor])
            reference_species_history = []
            time_history = []
            t = 0.0
            
            while t < estimated_IDT:  # Run until ignition delay is reached
                t = reactor_network.step()
                time_history.append(t)
                reference_species_history.append(gas[ref_species].X)
            i_ign = np.array(reference_species_history).argmax()  # Find index of maximum concentration of reference species
            tau = time_history[i_ign]  # Ignition delay time in seconds
            IDT[i][j] = tau * 1e6 # Convert to microseconds and store in selected array


# Calculate IDT for all equivalence ratios
print("Calculating IDT for phi = 0.5...")
calculate_Ignition_Delay(P_pa, T, gas, 0.5, store_arrays[0])

print("\nCalculating IDT for phi = 1.0...")
calculate_Ignition_Delay(P_pa, T, gas, 1.0, store_arrays[1])

print("\nCalculating IDT for phi = 2.0...")
calculate_Ignition_Delay(P_pa, T, gas, 2.0, store_arrays[2])

# -------------------- Plot Section Functions  --------------------------------

# Marker styles and colors to differentiate datasets
markers = ['o', 's', '^', 'D', 'v', '*', 'x', 'P']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:gray']


fig, axs = plt.subplots(3, 1, figsize=(10, 12))

def data_to_plot(exp_data, markers, colors, ax):
    num_columns = exp_data.shape[1]
    
    # Dynamically create labels based on column names
    pressure_labels = []
    for i in range(1, num_columns, 2):  # Get every second column name (the pressure columns)
        if i < num_columns:
            col_name = exp_data.columns[i]
            pressure_labels.append(col_name + " (Exp)")
    
    for i in range(0, num_columns, 2):
        if i + 1 < num_columns:  # Check if y column exists
            x = exp_data.iloc[:, i]
            y = exp_data.iloc[:, i + 1]
            
            # Remove NaN values and zeros (which cause plotting issues)
            mask = ~(pd.isna(x) | pd.isna(y) | (x == 0) | (y == 0))
            x_clean = x[mask]
            y_clean = y[mask]
            
            if len(x_clean) > 0:  # Only plot if there's data
                idx = i // 2
                marker = markers[idx % len(markers)]
                color = colors[idx % len(colors)]
                label = pressure_labels[idx] if idx < len(pressure_labels) else f"P{idx+1} (Exp)"
                ax.scatter(x_clean, y_clean, marker=marker, color=color, label=label, s=50) 

# Add secondary x-axis on the top for Temperature (K)
def T_to_invT(T):  # Converts T in K to 1000/T
    return 1000.0 / T

def invT_to_T(invT):  # Converts 1000/T back to T in K
    return 1000.0 / invT

# -------------------- Plot Section   --------------------------------
labels_pressure = ["1 atm", "3 atm", "5 atm", "10 atm", "25 atm", "44 atm"]
plt.tight_layout()
for plot_idx in range(3):
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))  
    
    # Plot experimental data
    data_to_plot(exp_datasets[plot_idx], markers, colors, ax)
    
    # Plot Cantera calculated data
    for i in range(len(P)):
        ax.semilogy(xaxis, store_arrays[plot_idx][i], 
                   linestyle='-', color=colors[i], 
                   label=labels_pressure[i] + ' (Cantera)')
    
    # Set labels, title, and formatting with increased font sizes
    ax.set_xlabel(r'$\frac{1000}{T}$ [1/K]', fontsize=16)
    ax.set_ylabel('Ignition Delay [μs]', fontsize=16)
    ax.set_title(f'Ignition Delay Time for φ = {phi[plot_idx]}', fontsize=18, pad=20)  
    ### ax.set_xlim(0.6, 1)  
    ### ax.set_ylim(bottom=1e2) 
    ax.legend(fontsize=12) 
    ax.grid(True)
    ax.tick_params(axis='both', which='major', labelsize=14)
    
 
    plt.tight_layout()
    # Add secondary x-axis
    secax = ax.secondary_xaxis('top', functions=(invT_to_T, T_to_invT))
    secax.set_xlabel('Temperature [K]', fontsize=16)
    secax.set_xticks(T)
    secax.tick_params(axis='x', labelsize=14)
    
    # Save individual figure
    filename = f"IDT_Phi_{phi[plot_idx]}_to_full.png"
    plt.savefig(os.path.join(save_path, filename), bbox_inches='tight', dpi=300)
    plt.close()