"""
This script was written for the purpose of assigments for the Combustion course on
Technical Univerity Munich.

Assigment: Excercise 4 - adabatic flame temp
Student: Herman Andersson

"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

equivalance_ratio = np.linspace(0.1, 2, 100)
T_store_AIR = [] # List to store adiabatic flame temperatures for each equivalence ratio
T_store_O2 = [] # List to store adiabatic flame temperatures for each equivalence ratio
gas = ct.Solution("gri30.yaml")
gas_O2 = ct.Solution("gri30.yaml")


for phi in equivalance_ratio:
    # Set the initial conditions and equivalence ratio for methane combustion
    gas.TP = 298, 101325 
    gas.set_equivalence_ratio(phi, 'CH4', 'O2:2, N2:7.52')
    
    gas_O2.TP = 298, 101325 
    gas_O2.set_equivalence_ratio(phi, 'CH4', 'O2:2')
    # Perform equilibrium calculation at constant enthalpy and pressure
    gas.equilibrate('HP', "auto")
    gas_O2.equilibrate('HP', "auto")
    
    # Store the adiabatic flame temperature for each equivalence ratio
    T_store_AIR.append(gas.T)
    T_store_O2.append(gas_O2.T)
    

#print(equivalance_ratio)
#print(T_store)

# Increase the font size of the plot
plt.rcParams.update({
    'axes.titlesize': 24,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'axes.labelsize': 28 
})

# Plot the adiabatic flame temperature against the equivalence ratio
plt.plot(equivalance_ratio, T_store_AIR)
plt.axhline(y=max(T_store_AIR), color='r', linestyle='--', label='Air Max Flame Temp')
plt.axvline(x=equivalance_ratio[np.argmax(T_store_AIR)], color='g', linestyle='--', label='Air phi')

plt.plot(equivalance_ratio, T_store_O2)
plt.axhline(y=max(T_store_O2), color='b', linestyle='--', label='O2 Max Flame Temp')
plt.axvline(x=equivalance_ratio[np.argmax(T_store_O2)], color='orange', linestyle='--', label='O2 phi')

plt.legend()
plt.xlabel('Equivalence Ratio (φ)')
plt.ylabel('Adiabatic Flame Temperature (K)')
plt.title('Adiabatic Flame Temperature vs Equivalence Ratio for Methane Combustion')

plt.grid()
plt.show()

print(f"Max Adiabatic Flame Temperature with air: {max(T_store_AIR):.2f} K at Equivalence Ratio: {equivalance_ratio[np.argmax(T_store_AIR)]:.2f},")
print(f"Max Adiabatic Flame Temperature with pure O2: {max(T_store_O2):.2f} K at Equivalence Ratio: {equivalance_ratio[np.argmax(T_store_O2)]:.2f},")