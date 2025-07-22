import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import re


def flame_sensitivity_analysis(pressure_atm, label, gri):
    """
    Perform flame speed sensitivity analysis for methane-air combustion
    """
    # Define gas phase mixture
    if gri == True:
        gas = ct.Solution("gri30.yaml")
    else:
        gas = ct.Solution("/home/brukare/Combustion/RMG_methane/cantera/chem.yaml")
       ### gas = ct.Solution("../cantera/chem.yaml") # Updated with temperature for range of 700-2000 K
        ### gas = ct.Solution("mechanisms_archive/rmg_99conv_1h25min_20250707/methane_99conv_1h25min.yaml")

    # Define initial conditions
    T0 = 298.0  # Initial temperature in K (298 K as per problem statement)
    P0 = pressure_atm * 1e5  # Convert bar to Pa
    
    # Set equivalence ratio for stoichiometric methane-air combustion
    if gri == True:
        gas.set_equivalence_ratio(1.1, 'CH4', 'O2:2.0, N2:7.52')
    else:
        gas.set_equivalence_ratio(1.1, 'CH4(1)', 'O2(2):2.0, N2:7.52')
    gas.TP = T0, P0
    
    # Domain width in meters
    width = 0.03
    
    # Create the flame object
    flame = ct.FreeFlame(gas, width=width)
    
    # Define tolerances for the solver
    flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    
    # Define logging level
    loglevel = 1
    
    # Solve for the flame
    flame.solve(loglevel=loglevel, auto=True)
    velocity = flame.velocity[0] * 100
    
    # Calculate flame speed sensitivities
    sensitivities = flame.get_flame_speed_reaction_sensitivities()
    
    # Get absolute values and sort
    abs_sensitivities = abs(sensitivities)
    sorted_indices = np.argsort(abs_sensitivities)[::-1]  # Sort in descending order
    
    # Filter reactions with sensitivity values >= 0.03
    reactions = []
    sensitivity_values = []
    
    for i in sorted_indices:
        if abs_sensitivities[i] >= 0.03:
            reactions.append(gas.reaction_equations()[i])
            sensitivity_values.append(sensitivities[i])
    
    return reactions, sensitivity_values, velocity

def normalize_reaction(reaction):
    """
    Remove species numbers and parentheses: H(5) + O2(2) -> H + O2
    """
    return re.sub(r'\([0-9]+\)', '', reaction).strip()

# Perform analysis for both pressure conditions
print("Performing sensitivity analysis at 1 bar for gri...")
reactions_gri, sens_gri, flame_speed_gri = flame_sensitivity_analysis(1.0, "Gri - 1.0 bar", True)

print("Performing sensitivity analysis at 1 bar for rmg...")
reactions_rmg, sens_rmg, flame_speed_rmg = flame_sensitivity_analysis(1.0, "Rmg - 1.0 bar", False)

# Create bar charts for reaction vs sensitivity value
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Calculate common x-axis limits
all_sens = sens_gri + sens_rmg
x_min = min(all_sens) * 1.1  # Add 10% margin
x_max = max(all_sens) * 1.1

# Plot for GRI - reverse order so largest sensitivity is at top
y_pos_1 = np.arange(len(reactions_gri))[::-1]  # Reverse order
colors_1 = ['red' if s > 0 else 'blue' for s in sens_gri]
bars1 = ax1.barh(y_pos_1, sens_gri, color=colors_1, alpha=0.7, height=0.8)
ax1.set_yticks(y_pos_1)
ax1.set_yticklabels(reactions_gri, fontsize=12)
ax1.set_xlabel('Sensitivity: ∂ln SL/∂ln k', fontsize=16)
ax1.set_ylabel('Reactions', fontsize=16)
ax1.set_title('GRI Mechanism\nFlame Speed Sensitivity Analysis\nφ = 1.1, p = 1.0 bar, T = 298 K', fontsize=16)
ax1.grid(True, axis='x', alpha=0.3)
ax1.axvline(x=0, color='black', linestyle='-', alpha=0.5)
ax1.tick_params(axis='x', labelsize=12)
ax1.set_xlim(x_min, x_max)

# Plot for RMG - reverse order so largest sensitivity is at top
y_pos_20 = np.arange(len(reactions_rmg))[::-1]  # Reverse order
colors_20 = ['red' if s > 0 else 'blue' for s in sens_rmg]
bars2 = ax2.barh(y_pos_20, sens_rmg, color=colors_20, alpha=0.7, height=0.8)
ax2.set_yticks(y_pos_20)
ax2.set_yticklabels(reactions_rmg, fontsize=12)
ax2.set_xlabel('Sensitivity: ∂ln SL/∂ln k', fontsize=16)
ax2.set_ylabel('Reactions', fontsize=16)
ax2.set_title('RMG Mechanism\nFlame Speed Sensitivity Analysis\nφ = 1.1, p = 1.0 bar, T = 298 K', fontsize=16)
ax2.grid(True, axis='x', alpha=0.3)
ax2.axvline(x=0, color='black', linestyle='-', alpha=0.5)
ax2.tick_params(axis='x', labelsize=12)
ax2.set_xlim(x_min, x_max)

plt.tight_layout()
plt.savefig('flame_speed_sensitivity_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# Print results
print(f"\nFlame speed for gri at 1 bar: {flame_speed_gri:.2f} cm/s")
print(f"Flame speed for rmg at 1 bar: {flame_speed_rmg:.2f} cm/s")

print(f"\nNumber of sensitive reactions for gri at 1 bar: {len(reactions_gri)}")
print(f"Number of sensitive reactions for rmg at 1 bar: {len(reactions_rmg)}")

# Print all sensitive reactions with their sensitivity values
print("\n" + "="*80)
print("ALL SENSITIVE REACTIONS FOR GRI AT 1 BAR")
print("="*80)
for i, (reaction, sensitivity) in enumerate(zip(reactions_gri, sens_gri)):
    print(f"{i+1:2d}. {reaction:<50} | Sensitivity: {sensitivity:+.4f}")

print("\n" + "="*80)
print("ALL SENSITIVE REACTIONS FOR RMG AT 1 BAR")
print("="*80)
for i, (reaction, sensitivity) in enumerate(zip(reactions_rmg, sens_rmg)):
    print(f"{i+1:2d}. {reaction:<50} | Sensitivity: {sensitivity:+.4f}")

# Simple comparison - normalize and compare
normalized_gri = [normalize_reaction(r) for r in reactions_gri]
normalized_rmg = [normalize_reaction(r) for r in reactions_rmg]

gri_set = set(normalized_gri)
rmg_set = set(normalized_rmg)

common = gri_set & rmg_set
only_gri = gri_set - rmg_set
only_rmg = rmg_set - gri_set

print(f"\nCommon reactions: {len(common)}")
print(f"Only in GRI: {len(only_gri)}")
print(f"Only in RMG: {len(only_rmg)}")

print("\n" + "="*60)
print("COMMON REACTIONS")
print("="*60)
for reaction in sorted(common):
    print(f"  {reaction}")

if only_gri:
    print("\n" + "="*60)
    print("ONLY IN GRI")
    print("="*60)
    for reaction in sorted(only_gri):
        print(f"  {reaction}")

if only_rmg:
    print("\n" + "="*60)
    print("ONLY IN RMG")
    print("="*60)
    for reaction in sorted(only_rmg):
        print(f"  {reaction}")