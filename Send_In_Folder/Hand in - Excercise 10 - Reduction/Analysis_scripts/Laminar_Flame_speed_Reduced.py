import pandas as pd
import numpy as np
import os

# needed to save the figures without displaying them
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

import matplotlib.pyplot as plt
import cantera as ct
# -------------------- Experimental Data import Section --------------------------------
# Import experimental flame speed data
script_dir = os.path.dirname(__file__)
os.chdir('/home/brukare/Combustion/RMG_methane')
exp_const_P = pd.read_csv(os.path.join(script_dir, 'FlameSpeedsExperimental_constant temperature_298K.csv'))
# Only take coloumns with relevant data
exp_phi = exp_const_P.iloc[:, 0] 
exp_LFS = exp_const_P.iloc[:, 1]    
# -------------------- Cantera Calculation Section --------------------------------
# Simulation Parameters
pressure = 1e5  # pressure of 1 bar in Pa
Tin_const = 298.0  # unburned gas temperature in K
width = 0.03  # width of the flame in m

# Solution objects for four mechanisms
gas_gri = ct.Solution("gri30.yaml")
gas_rmg = ct.Solution("/home/brukare/Combustion/RMG_methane/cantera/chem.yaml")
gas_R17 = ct.Solution("/home/brukare/Combustion/RMG_Reduced/reduced_17.yaml")
gas_R24 = ct.Solution("/home/brukare/Combustion/RMG_Reduced/reduced_24.yaml")

phi = np.linspace(0.5, 1.5, 21)  # equivalence ratio range
loglevel = 1  # amount of diagnostic output (0 to 8)

# Empty arrays to store results
LFS_gri = np.zeros(len(phi))  # GRI-Mech results
LFS_rmg= np.zeros(len(phi))  # RMG 99% (700-2000K)conv results
LFS_R17 = np.zeros(len(phi))  # GRI3 Reduced 17 results
LFS_R24 = np.zeros(len(phi))  # GRI3 Reduced 24 results

def calculate_flame_speed(phi, gas, loglevel, temperature, pressure, is_gri, result_array, width):
    """Calculate flame speed for a given mechanism"""
    for j in range(len(phi)):
        print("-"*50 + "*"*50 + "-"*50)
        print(f"Current iteration: {j}/{len(phi)-1}, phi = {phi[j]:.3f}")
        print("-"*50 + "*"*50 + "-"*50)
        
        # Set equivalence ratio based on mechanism type
        if is_gri:
            gas.set_equivalence_ratio(phi[j], 'CH4', 'O2:2, N2:7.52')
        else:
            gas.set_equivalence_ratio(phi[j], 'CH4(1)', 'O2(2):2, N2:7.52')
        
        # Set temperature and pressure
        gas.TP = temperature, pressure
        
        # Set up flame object
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=2.0, slope=0.06, curve=0.12)
        flame.solve(loglevel=loglevel, auto=True)

        # Store flame speed result in cm/s
        result_array[j] = flame.velocity[0] * 100

# Calculate flame speeds for all three mechanisms
print("Calculating GRI-Mech 3.0...")
calculate_flame_speed(phi, gas_gri, loglevel, Tin_const, pressure, True, LFS_gri, width)

print("\nCalculating for RMG 99% conv (700-2000K)...")
calculate_flame_speed(phi, gas_rmg, loglevel, Tin_const, pressure, False, LFS_rmg, width)

print("\nCalculating for Gri3 Reduced 17 ...")
calculate_flame_speed(phi, gas_R17, loglevel, Tin_const, pressure, True, LFS_R17, width)

print("\nCalculating for Gri3 Reduced 14 ...")
calculate_flame_speed(phi, gas_R24, loglevel, Tin_const, pressure, True, LFS_R24, width)
# -------------------- Plot Section --------------------------------

# Create comparison plot
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

ax.scatter(exp_phi, exp_LFS, marker='o', color='black', label='Experimental Data, 1atm', s=50)
# Plot all three mechanisms
ax.plot(phi, LFS_gri, linestyle='-', color='tab:blue', label='GRI-Mech 3.0')
ax.plot(phi, LFS_rmg, linestyle='-', color='tab:orange', label='RMG 99% conv (700-2000K)')
ax.plot(phi, LFS_R17, linestyle='-', color='tab:green', label='GRI3 Reduced 17')
ax.plot(phi, LFS_R24, linestyle='-', color='tab:red', label='GRI3 Reduced 24')

ax.set_title('Laminar Flame Speed Comparison at 298 K and 1 bar', fontsize=18)
ax.set_xlabel('Equivalence Ratio ϕ', fontsize=16)
ax.set_ylabel('Laminar Flame Speed [cm/s]', fontsize=16)
ax.legend(fontsize=14)
ax.grid(True)

plt.tight_layout()
plt.show()

# -------------------- Save Results --------------------------------

# Get current working directory
current_dir = os.getcwd()

# Create DataFrame with all results
results_df = pd.DataFrame({
    'phi': phi,
    'GRI_Mech_3.0': LFS_gri,
    'RMG 99% conv (700-2000K)': LFS_rmg,
    'GRI3 Reduced 17': LFS_R17,
    'GRI3 Reduced 24': LFS_R24
})

# Save to CSV in current directory
csv_filename = os.path.join(current_dir, 'flame_speed_comparison_reduced.csv')
results_df.to_csv(csv_filename, index=False)

# Save the figure in current directory
fig_filename = os.path.join(current_dir, 'FlameSpeed_Comparison_reduced.png')
plt.savefig(fig_filename, dpi=300, bbox_inches='tight')

print(f"Results saved to: {csv_filename}")
print(f"Figure saved to: {fig_filename}")
print("Calculation complete!")