import pandas as pd
import numpy as np
import time
import os

# needed to save the figures without displaying them
## import matplotlib
## matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

import matplotlib.pyplot as plt
import cantera as ct

# -------------------- Experimental Data import Section --------------------------------
# Import experimental flame speed data
script_dir = os.path.dirname(__file__)
os.chdir('/home/brukare/Combustion/RMG_methane')
IDT_1 = pd.read_csv(os.path.join(script_dir, 'IDT_1.0.csv'))

# -------------------- Cantera Calculation Section --------------------------------
# Simulation Parameters
phi = 1.1  # equivalence ratio
T_gri = np.linspace(700, 2000, 11)  # temperature in K // two present cause of troubleshooting, reamined incase of meachism errors
T_rmg = np.linspace(700, 2000, 11)  # temperature in K 
P_pa = 1e5  # pressure in pa (1 bar)
xaxis_gri = 1000/T_gri  # x-axis for the plot, 1000/T in K^-1
xaxis_rmg = 1000/T_rmg  # x-axis for the plot, 1000/T in K^-1
estimated_IDT = 1  # [s] Estimated ignition delay time in seconds 

# Define mechanism-specific species names, fuel and oxidiser
fuel_gri = 'CH4:1'     # GRI-Mech naming
fuel_rmg = 'CH4(1):1'  # RMG naming
ref_species_gri = 'CH'    # GRI-Mech naming  
ref_species_rmg = 'CH(9)' # RMG naming 

# Load mechanisms from backup folders
gas_gri = ct.Solution("gri30.yaml")
gas_rmg_fast = ct.Solution("mechanisms_archive/rmg_95conv_26min_20250707/methane_95conv_26min.yaml")
gas_rmg_slow = ct.Solution("mechanisms_archive/rmg_99conv_1h25min_20250707/methane_99conv_1h25min.yaml")
#gas_rmg = ct.Solution("cantera/chem.yaml")

ref_species = 'CH'  # Reference species for IDT determination

# Arrays to store results
idt_gri = np.zeros(len(T_gri))
idt_rmg_fast = np.zeros(len(T_rmg))
idt_rmg_slow = np.zeros(len(T_rmg))

def calculate_Ignition_Delay_simple(temperature_array, gas, mechanism_name, results_array, fuel):
    for j, temp in enumerate(temperature_array):
        print(f"Calculating {mechanism_name} - T={temp} K, P=1 bar, phi={phi}")
        
        # Create oxidiser string for phi
        O2_str = str(2 / phi)
        N2_str = str(7.52 / phi)
        if mechanism_name == "GRI-Mech 3.0":
            # For GRI-Mech, use the specific fuel and oxidiser format
            oxidiser = f'O2:{O2_str}, N2:{N2_str}'
        else:
            oxidiser = f'O2(2):{O2_str}, N2:{N2_str}'
        
        gas.set_equivalence_ratio(phi, fuel, oxidiser)
        gas.TP = temp, P_pa
        
        reactor = ct.Reactor(contents=gas)
        reactor_network = ct.ReactorNet([reactor])
        reference_species_history = []
        time_history = []
        t = 0.0
        
        while t < estimated_IDT:  # Run until ignition delay is reached
            t = reactor_network.step()
            time_history.append(t)
            if mechanism_name == "GRI-Mech 3.0":
                # For GRI-Mech, use the specific reference species
                reference_species_history.append(gas[ref_species_gri].X)
            else:
                # For RMG, use the specific reference species
                reference_species_history.append(gas[ref_species_rmg].X)
                
        i_ign = np.array(reference_species_history).argmax()  # Find index of maximum concentration of reference species
        tau = time_history[i_ign]  # Ignition delay time in seconds
        results_array[j] = tau * 1e6  # Convert to microseconds and store in selected array

# Calculate IDT for all three mechanisms
print("Calculating IDT for GRI-Mech 3.0...")
calculate_Ignition_Delay_simple(T_gri, gas_gri, "GRI-Mech 3.0", idt_gri, fuel_gri)

print("\nCalculating IDT for RMG 95% conv...")
calculate_Ignition_Delay_simple(T_rmg, gas_rmg_fast, "RMG 95% conv", idt_rmg_fast, fuel_rmg)

print("\nCalculating IDT for RMG 99% conv...")
calculate_Ignition_Delay_simple(T_rmg, gas_rmg_slow, "RMG 99% conv", idt_rmg_slow, fuel_rmg)
# -------------------- Plot Section --------------------------------
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
#  Plot experimental data
exp_temp = IDT_1.iloc[:, 0] 
exp_idt = IDT_1.iloc[:, 1]    
ax.scatter(exp_temp, exp_idt, marker='o', color='black', label='Experimental Data, 1atm, phi=1.0', s=50)
# Plot all three mechanisms
ax.semilogy(xaxis_gri, idt_gri, 'b-', linewidth=2, label='GRI-Mech 3.0')
ax.semilogy(xaxis_rmg, idt_rmg_fast, 'g-', linewidth=2,  label='RMG 95% conv (26min)')
ax.semilogy(xaxis_rmg, idt_rmg_slow, 'r-', linewidth=2, label='RMG 99% conv (1h25min)')

# Add secondary x-axis on the top for Temperature (K)
def invT_to_T(invT):  # Converts 1000/T back to T in K
    return 1000.0 / invT

def T_to_invT(T_val):  # Converts T in K to 1000/T
    return 1000.0 / T_val

# Set labels and formatting with increased font sizes
ax.set_xlabel(r'1000/T [K$^{-1}$]', fontsize=18)
ax.set_ylabel('Ignition Delay [μs]', fontsize=18)
ax.set_title(f'Ignition Delay Time Comparison - φ = {phi}, P = 1 bar', fontsize=20)
ax.legend(fontsize=14)
ax.set_xlim(0.5, 0.75) # Adjust x-axis limits for better visibility of experimental data  
ax.set_ylim(top=1e4) # Adjust y-axis limits for better visibility of experimental data  
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(True)

# Add secondary x-axis
secax = ax.secondary_xaxis('top', functions=(invT_to_T, T_to_invT))
secax.set_xlabel('Temperature [K]')

plt.tight_layout()
#plt.show()  # Show the plot in interactive mode
plt.savefig('/home/brukare/Combustion/RMG_methane/Analysis_scripts/ignition_delay_zoomed.png', dpi=300)

# Print summary
print("\n" + "="*60)
print("IGNITION DELAY COMPARISON SUMMARY")
print("="*60)
print(f"GRI-Mech 3.0:        {gas_gri.n_species} species, {gas_gri.n_reactions} reactions")
print(f"RMG 95% conv:        {gas_rmg_fast.n_species} species, {gas_rmg_fast.n_reactions} reactions")
print(f"RMG 99% conv:        {gas_rmg_slow.n_species} species, {gas_rmg_slow.n_reactions} reactions")
print("="*60)