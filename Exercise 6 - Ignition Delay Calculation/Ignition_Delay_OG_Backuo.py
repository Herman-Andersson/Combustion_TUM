import pandas as pd
import numpy as np
import time

# needed to save the figures without displaying them
#import matplotlib
#matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

import matplotlib.pyplot as plt
import cantera as ct

# -------------------- Experimental Data import Section --------------------------------
# Import experimental flame speed data
IDT_05 = pd.read_csv('/home/brukare/Combustion/Exercise 6 - Ignition Delay Calculation/IDT_0.5.csv')
IDT_1 = pd.read_csv('/home/brukare/Combustion/Exercise 6 - Ignition Delay Calculation/IDT_1.0.csv')
IDT_2 = pd.read_csv('/home/brukare/Combustion/Exercise 6 - Ignition Delay Calculation/IDT_2.0.csv')

# -------------------- Cantera Calculation Section --------------------------------
# Simulation Parameters
phi = [0.5,1.0,2.0]  # equivalence ratio range
T = np.linspace(1000, 2000, 10)  # temperature in K
P = [101325, 303975, 506625, 1013250]  # pressure in Pa
xaxis = 1000/T  # x-axis for the plot, 1000/T in K^-1

#######
estimated_IDT = 1 # see what value thisone should have 
#######

gas = ct.Solution("gri30.yaml")  # Solution object used to compute mixture properties
ref_species = 'CH'  # Reference species for IDT determination

# Empty arrays to store results
Store_IDT_05 = np.zeros((len(P), len(T)))  # Empty array to store
Store_IDT_1 = np.zeros((len(P), len(T)))  # Empty array to store
Store_IDT_2 = np.zeros((len(P), len(T)))  # Empty array to store

def calculate_Ignition_Delay(pressure, temperature, gas, phi, IDT):
    for i in range(len(pressure)):
        for j in range(len(temperature)):
            #print(f"Calculating IDT for P={pressure[i]/101325:.1f} atm, T={temperature[j]} K, phi={phi}")
            
            # Set equivalence ratio
            phi_02_str = "O2:" + str(2/phi)  # For phi = 0.5, we scale it down and convert to string
            phi_N2_str = "N2:" + str(7.52/phi)  # N2 is always in excess, so we scale it with phi
            oxidiser = phi_02_str + ', ' + phi_N2_str
            gas.set_equivalence_ratio(phi, 'CH4', oxidiser)
            gas.TP = temperature[j], pressure[i]  # Set temperature and pressure
            
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
            IDT[i][j] = tau


# Calculate IDT for all equivalence ratios
#print("Calculating IDT for phi = 0.5...")
calculate_Ignition_Delay(P, T, gas, 0.5, Store_IDT_05)

#print("\nCalculating IDT for phi = 1.0...")
calculate_Ignition_Delay(P, T, gas, 1.0, Store_IDT_1)

#print("\nCalculating IDT for phi = 2.0...")
calculate_Ignition_Delay(P, T, gas, 2.0, Store_IDT_2)

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
        ax.plot(x, y, linestyle='None', marker=marker, color=color, label=label)

# Add secondary x-axis on the top for Temperature (K)
def T_to_invT(T):  # Converts T in K to 1000/T
    return 1000.0 / T

def invT_to_T(invT):  # Converts 1000/T back to T in K
    return 1000.0 / invT

# -------------------- Plot Section   --------------------------------

# Plot first dataset (top subplot)
labels_presure = ["1 atm", "3 atm", "5 atm", "10 atm"]
data_to_plot(IDT_05, markers, colors, axs[0], labels_presure)
for i in range(len(P)):
    plots = axs[0].plot(xaxis[::-1], Store_IDT_05[i][::-1], marker=markers[i], linestyle='-', color=colors[i], label=labels_presure[i] + ' (Cantera)')
axs[0].set_xlabel(r'$\frac{1000}{T}$ [1/K]')
axs[0].set_ylabel('Ignition Delay [μs]')
axs[0].legend()
axs[0].grid(True)
secax = axs[0].secondary_xaxis('top', functions=(invT_to_T, T_to_invT))
secax.set_xlabel('Temperature [K]')
secax.set_xticks(T)

data_to_plot(IDT_1, markers, colors, axs[1], labels_presure)
axs[1].set_xlabel(r'$\frac{1000}{T}$ [1/K]')
axs[1].set_ylabel('Ignition Delay [μs]')
axs[1].legend()
axs[1].grid(True)
secax = axs[1].secondary_xaxis('top', functions=(invT_to_T, T_to_invT))
secax.set_xlabel('Temperature [K]')
secax.set_xticks(T)

data_to_plot(IDT_2, markers, colors, axs[2], labels_presure)
axs[2].set_xlabel(r'$\frac{1000}{T}$ [1/K]')
axs[2].set_ylabel('Ignition Delay [μs]')
axs[2].legend()
axs[2].grid(True)
secax = axs[2].secondary_xaxis('top', functions=(invT_to_T, T_to_invT))
secax.set_xlabel('Temperature [K]')
secax.set_xticks(T)


plt.tight_layout()
plt.show()

