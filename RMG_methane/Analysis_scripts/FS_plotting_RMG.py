import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

# -------------------- Experimental Data import Section --------------------------------
# Import experimental flame speed data
exp_const_temp = pd.read_csv('/home/brukare/Combustion/RMG_methane/Analysis_scripts/FlameSpeedsExperimental_constant temperature_298K.csv')
# Only take columns with relevant data
exp_phi = exp_const_temp.iloc[:, 0]  # Fixed: was exp_const_P
exp_LFS = exp_const_temp.iloc[:, 1]  # Fixed: was exp_const_P

# -------------------- Cantera Calculation import Section --------------------------------
# Import the saved flame speed results
flame_speed_results = pd.read_csv('/home/brukare/Combustion/RMG_methane/Analysis_scripts/flame_speed_comparison.csv')

# -------------------- Plot Section --------------------------------
# Create single plot as shown in the image
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Plot experimental data
ax.scatter(exp_phi, exp_LFS, marker='o', color='black', label='Experimental Data, 1atm', s=50)

# Plot the three mechanisms from the CSV file
ax.plot(flame_speed_results['phi'], flame_speed_results['GRI_Mech_3.0'], 
        linestyle='-', color='tab:blue', label='GRI-Mech 3.0')
ax.plot(flame_speed_results['phi'], flame_speed_results['RMG_95conv_26min'], 
        linestyle='-', color='tab:orange', label='RMG 95% conv (26 min)')
ax.plot(flame_speed_results['phi'], flame_speed_results['RMG_99conv_1h25min'], 
        linestyle='-', color='tab:green', label='RMG 99% conv (1h25min)')
ax.plot(flame_speed_results['phi'], flame_speed_results['RMG 99% conv (700-2000K)'], 
        linestyle='-', color='tab:red', label='RMG 99% conv (700-2000K)')

# Set plot properties
ax.set_title('Laminar Flame Speed Comparison at 298 K and 1 bar', fontsize=18)
ax.set_xlabel('Equivalence Ratio ϕ', fontsize=16)
ax.set_ylabel('Laminar Flame Speed [cm/s]', fontsize=16)
ax.legend(fontsize=14)
ax.grid(True)

plt.tight_layout()
plt.show()

# -------------------- Save Figure --------------------------------
# Get current working directory
import os
current_dir = os.getcwd()

# Save the figure in current directory
fig_filename = os.path.join(current_dir, 'FlameSpeed_Comparison_RMG.png')
plt.savefig(fig_filename, dpi=300, bbox_inches='tight')

print(f"Figure saved to: {fig_filename}")