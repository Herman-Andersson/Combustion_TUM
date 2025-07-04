# These values were taken from the previous assigment
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# Load AFT data 
AFT = np.loadtxt("Excercise 5 - NO formation pathways/AFT.txt", delimiter=',')
AFT = AFT.flatten()
Phi = np.linspace(0.1, 2.0, len(AFT))

# ------------------------------ Cantera calculation part ---------------------------------
# Code taken from Excercise 04, Problem 2
T_store_AIR = [] # List to store adiabatic flame temperatures for each equivalence ratio
gas = ct.Solution("gri30.yaml")

for e_ratio in Phi:
    # Set the initial conditions and equivalence ratio for methane combustion
    gas.TP = 298, 101325 
    gas.set_equivalence_ratio(e_ratio, 'CH4', 'O2:2, N2:7.52')
    # Perform equilibrium calculation at constant enthalpy and pressure
    gas.equilibrate('HP', "auto")
    # Store the adiabatic flame temperature for each equivalence ratio
    T_store_AIR.append(gas.T)

# ------------------------------ Cantera calculation end ---------------------------------

# Constants
R = 8.3145  # J/(mol·K)
P0 = 101325  # Pa
residence_time = 5e-3  # 5 ms in seconds
X_N2 = 0.79  # Initial mole fraction of N2
X_O2 = 0.21  # Initial mole fraction of O2
T_initial = 298  # K

def calculate_equilibrium_constant(T):
    # Calculate equilibrium constant for O2 ⇌ 2O reaction
    theta = T/1000
    a1, a2, a3, a4, a5, a6 = 6.434, -0.2755, 0.02396, -0.111e-2, 0.8258, -25.80
    log10_Kp = a1 + a2*theta + a3*theta**2 + a4*theta**3 + a5*np.log(theta) + a6/theta
    Kp = 10**log10_Kp
    return Kp

def calculate_rate_constants(T):
    # forward rate constants for Zeldovich mechanism
    k2f = 1.8e14 * np.exp(-318000/(R*T)) * 1e-6 # m3/(mol·s), cm3 = m3 * 1e-6
    k3f = 9e9 * np.exp(-27000/(R*T)) * 1e-6 # m3/(mol·s), cm3 = m3 * 1e-6
    return k2f, k3f

def molar_frac_to_concentration(X, T):
	# Convert mole fraction to molar concentration
	# Species concentration (in ppm) = x_species * 10**6
	return X * P0 / (R * T)  # mol/m³

def concentration_to_molar_frac(concenctration, T):
	# Convert molar concentration to mole fraction
	# Species concentration (in ppm) = x_species * 10**6
	return concenctration*R*T / P0

# Might not be needed
#def O_concentration(T):
#    Kp = calculate_equilibrium_constant(T)
#    O2_conc = molar_frac_to_concentration(X_O2, T)  # O2 concentration in mol/m³
#    O_conc = (O2_conc *((Kp*P0)/(R*T)**0.5))  # O concentration in mol/m³
#    return O_conc

def NO_concentration(T):
    N2_conc = molar_frac_to_concentration(X_N2, T)  # N2 concentration in mol/m³
    O2_conc = molar_frac_to_concentration(X_O2, T)  # N2 concentration in mol/m³
    Kp = calculate_equilibrium_constant(T)  
    k2f, k3f = calculate_rate_constants(T)  # Get rate constants

    dNO_dt = 2 * k2f * ((Kp*P0)/(R*T))**0.5 * N2_conc * (O2_conc)**0.5

    NO_conc = dNO_dt * residence_time  # mol/m³
    return NO_conc

def convert_to_ppm(fraction, T):
	# Convert concentration in mol/m³ to ppm
	ppm_conc = fraction * 1e6
	return ppm_conc

NO_conc_manual = np.zeros(len(Phi))
NO_frac_manual = np.zeros(len(Phi))
NO_ppm_manual = np.zeros(len(Phi))

NO_conc_cantera = np.zeros(len(Phi))
NO_frac_cantera = np.zeros(len(Phi))
NO_ppm_cantera = np.zeros(len(Phi))

for i in range(len(Phi)):
    NO_conc_manual[i] = NO_concentration(AFT[i])
    NO_frac_manual[i] = concentration_to_molar_frac(NO_conc_manual[i], AFT[i])
    NO_ppm_manual[i] = convert_to_ppm(NO_frac_manual[i], AFT[i])

    NO_conc_cantera[i] = NO_concentration(T_store_AIR[i])
    NO_frac_cantera[i] = concentration_to_molar_frac(NO_conc_cantera[i], T_store_AIR[i])
    NO_ppm_cantera[i] = convert_to_ppm(NO_frac_cantera[i], T_store_AIR[i])


# Print part to see function of code above works
#print("NO Concentration (mol/m³):", NO_conc)
#print("NO Concentration (ppm):", NO_ppm)
#print("AFT (K):", AFT)

# Left plot with Manual AFT and NO concentration
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16,6))
Color_B = "tab:blue"
Color_R = "tab:red"

ax1.set_xlabel('Equivalence ratio, φ', fontsize=12)
ax1.set_ylabel('Adiabatic flame temperature [K]', color=Color_R, fontsize=12)
ax1.plot(Phi, AFT, color=Color_R, label='AFT')
ax1.axhline(y=max(AFT), color="orange", linestyle='--', label='AFT Max')
ax1.tick_params(axis='y', labelcolor=Color_R)

ax1b = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax1b.set_ylabel('NO Concentration (ppm)', color=Color_B, fontsize=12)
ax1b.plot(Phi, NO_ppm_manual, color=Color_B, label='NO (ppm)')
ax1b.axhline(y=max(NO_ppm_manual), color="black", linestyle='--', label='NO (ppm) Max')
ax1b.tick_params(axis='y', labelcolor=Color_B)

# ----------------------------------------------------------------
# ----------------------------------------------------------------

# Right plot with Cantera AFT and NO concentration
ax2.set_xlabel('Equivalence ratio, φ', fontsize=12)
ax2.set_ylabel('Adiabatic flame temperature [K] Cantera', color=Color_R, fontsize=12)
ax2.plot(Phi, T_store_AIR, color=Color_R, label='AFT')
ax2.axhline(y=max(T_store_AIR), color="orange", linestyle='--', label='AFT Max')
ax2.tick_params(axis='y', labelcolor=Color_R)

ax2b = ax2.twinx()  # instantiate a second axes that shares the same x-axis

ax2b.set_ylabel('NO Concentration (ppm)', color=Color_B, fontsize=12)
ax2b.plot(Phi, NO_ppm_cantera, color=Color_B, label='NO (ppm)')
ax2b.axhline(y=max(NO_ppm_cantera), color="black", linestyle='--', label='NO (ppm) Max')
ax2b.tick_params(axis='y', labelcolor=Color_B)

fig.suptitle('Methane Combustion: AFT and NO Formation', fontsize=16)
fig.tight_layout()
plt.show()

print("-" * 10 + "\ncalculation completed\n" + "-" * 10)

# Legends?

# -----------------------------------------------
# Test to see if volumes are the same, remove later
# Calculate the volume for the AFT and Cantera solutions (assuming 1 mol of gas)
# V = nRT/P, with n=1 mol

volume_AFT = R * AFT / P0  # m³ for each AFT temperature
volume_cantera = R * np.array(T_store_AIR) / P0  # m³ for each Cantera temperature

# Find index where NO ppm is at the highest value for both cases
idx_max_NO_manual = np.argmax(NO_ppm_manual)
idx_max_NO_cantera = np.argmax(NO_ppm_cantera)

print("Volume for AFT solution (m³):", volume_AFT)
print("Volume for Cantera solution (m³):", volume_cantera)

print(f"\nAt max NO (manual):")
print(f"  Phi = {Phi[idx_max_NO_manual]:.3f}")
print(f"  Volume = {volume_AFT[idx_max_NO_manual]:.6f} m³")

print(f"\nAt max NO (cantera):")
print(f"  Phi = {Phi[idx_max_NO_cantera]:.3f}")
print(f"  Volume = {volume_cantera[idx_max_NO_cantera]:.6f} m³")

if np.isclose(volume_AFT[idx_max_NO_manual], volume_cantera[idx_max_NO_cantera], rtol=1e-6):
    print("\nThe volumes at maximum NO (ppm) are the same (within tolerance).")
else:
    print("\nThe volumes at maximum NO (ppm) are different.")