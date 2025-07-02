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
IDT_05 = pd.read_csv(os.path.join(script_dir, 'IDT_0.5.csv'))
IDT_1 = pd.read_csv(os.path.join(script_dir, 'IDT_1.0.csv'))
IDT_2 = pd.read_csv(os.path.join(script_dir, 'IDT_2.0.csv'))
exp_datasets = [IDT_05, IDT_1, IDT_2]
# -------------------- Cantera Calculation Section --------------------------------
# Simulation Parameters
phi = [0.5,1.0,2.0]  # equivalence ratio range
T = np.linspace(1000, 2000, 10)  # temperature in K
P = [1, 3, 5, 10]  # pressure in atm
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

def data_to_plot(exp_data, markers, colors, ax, labels):
    num_columns = exp_data.shape[1]
    for i in range(0, num_columns, 2):
        x = exp_data.iloc[:, i]
        y = exp_data.iloc[:, i + 1]
        idx = i // 2
        marker = markers[idx % len(markers)]
        color = colors[idx % len(colors)]
        label = labels[idx % len(labels)] + " (Exp)" 
        ax.scatter(x, y, marker=marker, color=color, label=label, s=50)  # BONUS: Use scatter instead of plot for exp data

# Add secondary x-axis on the top for Temperature (K)
def T_to_invT(T):  # Converts T in K to 1000/T
    return 1000.0 / T

def invT_to_T(invT):  # Converts 1000/T back to T in K
    return 1000.0 / invT

# -------------------- Plot Section   --------------------------------
labels_pressure = ["1 atm", "3 atm", "5 atm", "10 atm"]

for plot_idx in range(3):
    # Plot experimental data
    data_to_plot(exp_datasets[plot_idx], markers, colors, axs[plot_idx], labels_pressure)
    
    # Plot Cantera calculated data
    for i in range(len(P)):
        axs[plot_idx].semilogy(xaxis, store_arrays[plot_idx][i], 
                              linestyle='-', color=colors[i], 
                              label=labels_pressure[i] + ' (Cantera)')
    
    # Set labels, title, and formatting
    axs[plot_idx].set_xlabel(r'$\frac{1000}{T}$ [1/K]')
    axs[plot_idx].set_ylabel('Ignition Delay [μs]')
    axs[plot_idx].set_title(f'Ignition Delay Time for φ = {phi[plot_idx]}')
    axs[plot_idx].legend()
    axs[plot_idx].grid(True)
    
    # Add secondary x-axis
    secax = axs[plot_idx].secondary_xaxis('top', functions=(invT_to_T, T_to_invT))
    secax.set_xlabel('Temperature [K]')
    secax.set_xticks(T)

plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), "IDT_Combined_Various_phi.png"))
#plt.show()