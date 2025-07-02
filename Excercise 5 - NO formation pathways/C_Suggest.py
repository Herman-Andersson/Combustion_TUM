# NO concentration prediction using Zeldovich mechanism
# Exercise 5: Prediction of NO concentration at varying equivalence ratio

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Given data from previous assignment
Phi = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
AFT = np.array([1550.43066213, 1707.93580174, 1857.09049583, 1998.88048322, 2134.07544956,
                2033.27743736, 1973.03478179, 1899.11353857, 1828.80502939, 1761.62669197,
                1697.19220524, 1635.19437898, 1575.38622149, 1517.5654539, 1461.56302568])

# Constants
R = 8.3145  # J/(mol·K)
P0 = 101325  # Pa
residence_time = 5e-3  # 5 ms in seconds
X_N2 = 0.79  # Initial mole fraction of N2
X_O2 = 0.21  # Initial mole fraction of O2
T_initial = 298  # K

def calculate_equilibrium_constant(T):
    """Calculate equilibrium constant for O2 ⇌ 2O reaction"""
    theta = T / 1000
    a1, a2, a3, a4, a5, a6 = 6.434, -0.2755, 0.02396, -0.111e-2, 0.8258, -25.80
    
    log10_Kp = a1 + a2*theta + a3*theta**2 + a4*theta**3 + a5*np.log(theta) + a6/theta
    Kp = 10**log10_Kp
    return Kp

def calculate_rate_constants(T):
    """Calculate forward rate constants for Zeldovich mechanism"""
    # Convert from cm³/(mol·s) to m³/(mol·s) by multiplying by 1e-6
    k2f = 1.8e14 * np.exp(-318000/(R*T)) * 1e-6  # m³/(mol·s)
    k3f = 9.0e9 * np.exp(-27000/(R*T)) * 1e-6    # m³/(mol·s)
    return k2f, k3f

def calculate_species_concentrations(phi, T):
    """Calculate molar concentrations of species at equilibrium"""
    # For methane combustion: CH4 + 2(O2 + 3.76N2) → CO2 + 2H2O + 7.52N2
    # At stoichiometric: 1 mol CH4 + 2 mol O2 + 7.52 mol N2
    
    # Calculate mole fractions after combustion (simplified approach)
    # Assuming complete combustion and using initial air composition
    
    # Molar concentration conversion: [species] = X_species * P / (R * T)
    total_pressure = P0  # Assuming constant pressure
    
    # For simplification, using initial air composition scaled by temperature
    # In real combustion, these would change based on equilibrium calculations
    conc_N2 = X_N2 * total_pressure / (R * T)  # mol/m³
    conc_O2 = X_O2 * total_pressure / (R * T)  # mol/m³
    
    return conc_N2, conc_O2

def calculate_NO_concentration(phi, T):
    """Calculate NO concentration using Zeldovich mechanism"""
    if T < 1000:  # NO formation negligible at low temperatures
        return 0.0
    
    # Get rate constants
    k2f, k3f = calculate_rate_constants(T)
    
    # Get equilibrium constant
    Kp = calculate_equilibrium_constant(T)
    
    # Get species concentrations
    conc_N2, conc_O2 = calculate_species_concentrations(phi, T)
    
    # Calculate [O] from equilibrium: O2 ⇌ 2O
    # Kp = [O]² / [O2] * (P0 / (R*T))
    # [O] = sqrt(Kp * [O2] * P0 / (R*T))
    conc_O = np.sqrt(Kp * conc_O2 * P0 / (R * T))
    
    # Calculate NO formation rate using equation (11)
    # d[NO]/dt = 2 * k2f * sqrt(Kp * P0 / (R*T)) * [N2] * [O2]^(1/2)
    dNO_dt = 2 * k2f * np.sqrt(Kp * P0 / (R * T)) * conc_N2 * np.sqrt(conc_O2)
    
    # Calculate NO concentration after residence time
    NO_conc = dNO_dt * residence_time  # mol/m³
    
    return NO_conc

def convert_to_ppm(NO_conc, T):
    """Convert NO concentration from mol/m³ to ppm"""
    # Total molar concentration at temperature T
    total_conc = P0 / (R * T)  # mol/m³
    
    # Mole fraction
    X_NO = NO_conc / total_conc
    
    # Convert to ppm
    NO_ppm = X_NO * 1e6
    
    return NO_ppm

# Calculate NO concentrations for all equivalence ratios
NO_concentrations = []
NO_ppm_values = []

for i, (phi, T) in enumerate(zip(Phi, AFT)):
    NO_conc = calculate_NO_concentration(phi, T)
    NO_ppm = convert_to_ppm(NO_conc, T)
    
    NO_concentrations.append(NO_conc)
    NO_ppm_values.append(NO_ppm)
    
    print(f"φ = {phi:.1f}, T = {T:.1f} K, NO = {NO_ppm:.2f} ppm")

NO_concentrations = np.array(NO_concentrations)
NO_ppm_values = np.array(NO_ppm_values)

# Create plots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

# Plot 1: Adiabatic Flame Temperature
ax1.plot(Phi, AFT, 'b-o', linewidth=2, markersize=6, label='AFT')
ax1.axhline(y=max(AFT), color='r', linestyle='--', alpha=0.7, label=f'Max AFT = {max(AFT):.0f} K')
ax1.axhline(y=1850, color='orange', linestyle=':', alpha=0.7, label='1850 K (NO formation threshold)')
ax1.set_xlabel('Equivalence ratio, φ')
ax1.set_ylabel('Adiabatic flame temperature [K]')
ax1.set_title('Adiabatic Flame Temperature vs Equivalence Ratio')
ax1.grid(True, alpha=0.3)
ax1.legend()

# Plot 2: NO concentration in ppm
ax2.plot(Phi, NO_ppm_values, 'r-o', linewidth=2, markersize=6, label='NO concentration')
ax2.set_xlabel('Equivalence ratio, φ')
ax2.set_ylabel('NO concentration [ppm]')
ax2.set_title('NO Concentration vs Equivalence Ratio')
ax2.grid(True, alpha=0.3)
ax2.legend()

# Plot 3: Combined plot showing both AFT and NO
ax3_temp = ax3.twinx()

line1 = ax3.plot(Phi, NO_ppm_values, 'r-o', linewidth=2, markersize=6, label='NO concentration')
line2 = ax3_temp.plot(Phi, AFT, 'b-s', linewidth=2, markersize=6, label='AFT')

ax3.set_xlabel('Equivalence ratio, φ')
ax3.set_ylabel('NO concentration [ppm]', color='r')
ax3_temp.set_ylabel('Adiabatic flame temperature [K]', color='b')
ax3.set_title('NO Concentration and AFT vs Equivalence Ratio')
ax3.grid(True, alpha=0.3)

# Combine legends
lines1, labels1 = ax3.get_legend_handles_labels()
lines2, labels2 = ax3_temp.get_legend_handles_labels()
ax3.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.tight_layout()
plt.show()

# Analysis and comments
print("\n" + "="*60)
print("ANALYSIS OF RESULTS")
print("="*60)

print(f"\nTemperature range for significant NO formation:")
print(f"- Threshold temperature: ~1850 K")
print(f"- Maximum AFT: {max(AFT):.1f} K at φ = {Phi[np.argmax(AFT)]:.1f}")
print(f"- Maximum NO: {max(NO_ppm_values):.2f} ppm at φ = {Phi[np.argmax(NO_ppm_values)]:.1f}")

print(f"\nKey observations:")
print(f"- NO formation is negligible below 1850 K")
print(f"- Peak NO formation occurs near stoichiometric conditions")
print(f"- NO decreases with increasing φ due to lower oxygen availability")
print(f"- NO also decreases with decreasing φ due to lower temperatures")

# Find equivalence ratios where AFT > 1850 K
high_temp_indices = AFT > 1850
print(f"\nEquivalence ratios with AFT > 1850 K:")
for phi, aft, no_ppm in zip(Phi[high_temp_indices], AFT[high_temp_indices], NO_ppm_values[high_temp_indices]):
    print(f"  φ = {phi:.1f}: AFT = {aft:.1f} K, NO = {no_ppm:.2f} ppm")