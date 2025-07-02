import cantera as ct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# Chemical reaction
fuel = 'CH4:1'
mech = 'gri30.yaml'  # Mechanism file

# Variables
pressures = [1, 3, 5, 10]  # atm
phis = [0.5] # [0.5, 1.0, 2.0]  # Equivalence ratios
temperatures = np.linspace(1000, 2000, 10)  # in K
xaxis = 1000 / temperatures

# Reference species to check ignition
ref_species = 'CH'

# Convert pressures to Pascals
pressures_pa = [p * ct.one_atm for p in pressures]

# Function to calculate ignition delay time
def ignition_delay(gas, T, P):
    gas.TP = T, P
    reactor = ct.Reactor(contents=gas)
    reactor_network = ct.ReactorNet([reactor])
    time = 0.0
    estimated_IDT = 1.0  # seconds
    reference_species_history = []
    time_history = []
    while time < estimated_IDT:
        time = reactor_network.step() # Advance the time
        time_history.append(time)
        reference_species_history.append(gas[ref_species].X)
    i_ign = np.array(reference_species_history).argmax() # Find index at which ref_species is highest
    return time_history[i_ign] # return ignition delay time

# Function to read the data from the excel file
def read_data(i, phi):
    filename = "IDT_" + str(phi) + ".csv"
    print(filename)
    P = 'PHI' + str(i+1)
    atm = 'IDT_' + str(pressures[i]) + 'atm(microsec)'
    x_axis = pd.read_csv(filename, usecols=[P])
    y_axis = pd.read_csv(filename, usecols=[atm])
    return [x_axis, y_axis]

# Create location to store the results
IDT = {}

for phi in phis:
    IDT[phi] = {}
    for P, P_atm in zip(pressures_pa, pressures):
        delay_times = []
        for T in temperatures:
            gas = ct.Solution(mech)
            # Create oxidiser string for phi
            O2 = 'O2:' + str(2 / phi)
            N2 = 'N2:' + str(7.52 / phi)
            oxidiser = O2 + ', ' + N2
            gas.set_equivalence_ratio(phi, fuel, oxidiser)
            delay = ignition_delay(gas, T, P) * 1e6 # micro s
            delay_times.append(delay)
            print(f"phi={phi}, P={P_atm} atm, T={T} K -> τ = {delay:.5f} mu s")
        IDT[phi][P_atm] = delay_times


fig, ax1 = plt.subplots(figsize=(8, 5))

colors = ['k', 'r', 'b', 'g']
linestyles = ['-', '--', '-.', ':']
markers = ['.', 'x', 'D', '^']

for i, phi in enumerate(phis):
    for j, P in enumerate(pressures):
        data = read_data(j, phi)
        label_cantera = f'{P} atm (Cantera)'
        label_experimentaly = f'{P} atm (Exp)'
        taus = IDT[phi][P]
        ax1.semilogy(xaxis, taus, color=colors[j % 4], label=label_cantera)
        ax1.scatter(data[0], data[1], color=colors[j % 4], marker=markers[j % 4], label=label_experimentaly)

ax1.set_xlabel('1000 / T [1/K]')
ax1.set_ylabel('Ignition Delay Time ' + r'[$\mu$s]')
ax1.grid(True)
ax1.legend()

# Secondary x-axis for T
def x1_to_x2(x): return 1000 / x

def x2_to_x1(x): return 1000 / x

ax2 = ax1.secondary_xaxis('top', functions=(x1_to_x2, x2_to_x1))
ax2.set_xlabel('Temperature [K]')

#plt.title('Ignition Delay Times for Various Conditions')
plt.tight_layout()
plt.savefig('Exercise_6.png')
