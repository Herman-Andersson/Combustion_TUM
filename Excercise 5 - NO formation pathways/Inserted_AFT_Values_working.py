# These values were taken from the previous assigment
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
 
# Might change to actual calculated values in code, this for testing the code in generall (04/06/2025)
Phi = np.arange(0.6, 2.1, 0.1)
AFT = [1550.43066213, 1707.93580174, 1857.09049583, 1998.88048322, 2134.07544956,
	   2033.27743736, 1973.03478179, 1899.11353857, 1828.80502939, 1761.62669197,
	   1697.19220524, 1635.19437898, 1575.38622149, 1517.5654539, 1461.56302568]


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

NO_conc = np.zeros(len(Phi))
NO_frac = np.zeros(len(Phi))
NO_ppm = np.zeros(len(Phi))

for i in range(len(Phi)):
    NO_conc[i] = NO_concentration(AFT[i])
    NO_frac[i] = concentration_to_molar_frac(NO_conc[i], AFT[i])
    NO_ppm[i] = convert_to_ppm(NO_frac[i], AFT[i])

# Used to see code values
#print("NO Concentration (mol/m³):", NO_conc)
#print("NO Concentration (ppm):", NO_ppm)
#print("AFT (K):", AFT)

# Plotting AFT and NO_ppm against Phi in one plot with twin y-axes
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Equivalence ratio, φ', fontsize=12)
ax1.set_ylabel('Adiabatic flame temperature [K]', color=color, fontsize=12)
ax1.plot(Phi, AFT, color=color, marker='o', label='AFT')
ax1.axhline(y=max(AFT), color="orange", linestyle='--', label='AFT Max')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('NO Concentration (ppm)', color=color, fontsize=12)
ax2.plot(Phi, NO_ppm, color=color, marker='s', label='NO (ppm)')
ax2.axhline(y=max(NO_ppm), color="black", linestyle='--', label='NO (ppm) Max')
ax2.tick_params(axis='y', labelcolor=color)

fig.suptitle('Methane Combustion: AFT and NO Formation', fontsize=16)
fig.tight_layout()
plt.show()

print("Calculation completed!")
