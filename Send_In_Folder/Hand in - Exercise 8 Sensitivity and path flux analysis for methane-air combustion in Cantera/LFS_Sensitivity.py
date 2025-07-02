import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

def flame_sensitivity_analysis(pressure_atm, label):
    """
    Perform flame speed sensitivity analysis for methane-air combustion
    """
    # Define gas phase mixture
    gas = ct.Solution("gri30.yaml")
    
    # Define initial conditions
    T0 = 298.0  # Initial temperature in K (298 K as per problem statement)
    P0 = pressure_atm * ct.one_atm  # Convert atm to Pa
    
    # Set equivalence ratio for stoichiometric methane-air combustion
    gas.set_equivalence_ratio(1.0, 'CH4', 'O2:2.0, N2:7.52')
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

# Perform analysis for both pressure conditions
print("Performing sensitivity analysis at 1 atm...")
reactions_1atm, sens_1atm, flame_speed_1atm = flame_sensitivity_analysis(1.0, "1.0 atm")

print("Performing sensitivity analysis at 20 atm...")
reactions_20atm, sens_20atm, flame_speed_20atm = flame_sensitivity_analysis(20.0, "20.0 atm")

# Create bar charts for reaction vs sensitivity value
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Plot for 1 atm - reverse order so largest sensitivity is at top
y_pos_1 = np.arange(len(reactions_1atm))[::-1]  # Reverse order
colors_1 = ['red' if s > 0 else 'blue' for s in sens_1atm]
bars1 = ax1.barh(y_pos_1, sens_1atm, color=colors_1, alpha=0.7, height=0.8)
ax1.set_yticks(y_pos_1)
ax1.set_yticklabels(reactions_1atm, fontsize=12)
ax1.set_xlabel('Sensitivity: ∂ln SL/∂ln k', fontsize=14)
ax1.set_ylabel('Reactions', fontsize=14)
ax1.set_title('Flame Speed Sensitivity Analysis\nφ = 1.0\np = 1.0 atm\nT = 298 K', fontsize=12)
ax1.grid(True, axis='x', alpha=0.3)
ax1.axvline(x=0, color='black', linestyle='-', alpha=0.5)
ax1.tick_params(axis='x', labelsize=12)

# Plot for 20 atm - reverse order so largest sensitivity is at top
y_pos_20 = np.arange(len(reactions_20atm))[::-1]  # Reverse order
colors_20 = ['red' if s > 0 else 'blue' for s in sens_20atm]
bars2 = ax2.barh(y_pos_20, sens_20atm, color=colors_20, alpha=0.7, height=0.8)
ax2.set_yticks(y_pos_20)
ax2.set_yticklabels(reactions_20atm, fontsize=12)
ax2.set_xlabel('Sensitivity: ∂ln SL/∂ln k', fontsize=14)
ax2.set_ylabel('Reactions', fontsize=14)
ax2.set_title('Flame Speed Sensitivity Analysis\nφ = 1.0\np = 20.0 atm\nT = 298 K', fontsize=12)
ax2.grid(True, axis='x', alpha=0.3)
ax2.axvline(x=0, color='black', linestyle='-', alpha=0.5)
ax2.tick_params(axis='x', labelsize=12)

plt.tight_layout()
plt.savefig('flame_speed_sensitivity_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# Print results
print(f"\nFlame speed at 1 atm: {flame_speed_1atm:.2f} cm/s")
print(f"Flame speed at 20 atm: {flame_speed_20atm:.2f} cm/s")

print(f"\nNumber of sensitive reactions at 1 atm: {len(reactions_1atm)}")
print(f"Number of sensitive reactions at 20 atm: {len(reactions_20atm)}")

# Print all sensitive reactions with their sensitivity values
print("\n" + "="*80)
print("ALL SENSITIVE REACTIONS AT 1 ATM")
print("="*80)
for i, (reaction, sensitivity) in enumerate(zip(reactions_1atm, sens_1atm)):
    print(f"{i+1:2d}. {reaction:<50} | Sensitivity: {sensitivity:+.4f}")

print("\n" + "="*80)
print("ALL SENSITIVE REACTIONS AT 20 ATM")
print("="*80)
for i, (reaction, sensitivity) in enumerate(zip(reactions_20atm, sens_20atm)):
    print(f"{i+1:2d}. {reaction:<50} | Sensitivity: {sensitivity:+.4f}")

# Find common and different reactions
reactions_1atm_set = set(reactions_1atm)
reactions_20atm_set = set(reactions_20atm)

common_reactions = reactions_1atm_set.intersection(reactions_20atm_set)
unique_1atm = reactions_1atm_set - reactions_20atm_set
unique_20atm = reactions_20atm_set - reactions_1atm_set

print(f"\nCommon sensitive reactions: {len(common_reactions)}")
print(f"Reactions only sensitive at 1 atm: {len(unique_1atm)}")
print(f"Reactions only sensitive at 20 atm: {len(unique_20atm)}")

# Analysis of differences
print("\n" + "="*60)
print("ANALYSIS OF PRESSURE EFFECTS ON SENSITIVITY")
print("="*60)

if len(unique_1atm) > 0:
    print("\nReactions only sensitive at 1 atm:")
    for reaction in unique_1atm:
        print(f"  - {reaction}")

if len(unique_20atm) > 0:
    print("\nReactions only sensitive at 20 atm:")
    for reaction in unique_20atm:
        print(f"  - {reaction}")

