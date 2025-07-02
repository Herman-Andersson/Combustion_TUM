import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define your fuel, methane (CH4) coefficients
x = 1   #insert amount of carbon in fuel
y = 4   #insert amount of hydrogen in fuel

equivalance_ratio = np.linspace(0.1, 2, 100)

def calculate_CP(T,coeffs):
    """Calculate specific heat capacity using polynomial coefficients"""
    a1, a2, a3, a4, a5 = coeffs[:5]
    Ru = 8.314  # J/(mol·K)
    cp_bar = Ru * (a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4)
    return cp_bar / 1000  # Convert to kJ/(kmol·K)

def water_gas_equilibrium_kp(T):
    """Calculate equilibrium constant for water-gas shift reaction"""
    return 10**(-2.4198 + 0.0003855*T + 2180.6/T)

def solve_species_balance(phi, T):
    """Solve species balance for given equivalence ratio and temperature"""
    a = (x + y/4) / phi  # Stoichiometric oxygen requirement    
    Kp = water_gas_equilibrium_kp(T)
    if phi <= 1.0: # Lean mixture
        # CxHy + a*O2 + 3.76*a*N2 -> b*CO2 + d*H2O + f*O2 + 3.76*a*N2
        b = x
        d = y/2
        f = a - b - d/2
        c = 0 # No CO
        e = 0 # No H2, all to H2O
    else: # Rich mixture
        # CxHy + a*O2 + 3.76*a*N2 -> b*CO2 + c*CO + d*H2O + e*H2 + f*O2 + 3.76*a*N2 
        # c = x - b  (equation 1)
        # d = 2*a - b - x  (equation 2)
        # e = -2*a + b + x + y/2  (equation 3)
        # Kp = (b*e)/(c*d)  (equation 4)

        def equilibrium_eq(b): # probably need antoher method cause the guess determines everything
            c = x - b                      # CO
            d = 2*a - b - x               # H2O
            e = -2*a + b + x + y/2        # H2
            if c <= 0 or d <= 0:          # avoid division by zero or negative mols
                return 1e6  # Penalize invalid solutions
            return (b * e) / (c * d) - Kp

        # Initial guess for b
        b_guess = x * 0.5
        b = fsolve(equilibrium_eq, b_guess)[0]

        # Back-calculate remaining species
        c = x - b
        d = 2*a - b - x
        e = -2*a + b + x + y/2
        f = 0  # Assume no excess O2 in rich case

    return {'a': a, 'b': b, 'c': c, 'd': d, 'e': e, 'f': f, 'Kp': Kp}

#--------------------
# Example test
phi_test = 1.3
T_test = 1000  # in Kelvin
result = solve_species_balance(phi_test, T_test)
print(f"phi = {phi_test}, T = {T_test} K")
print(result)
#--------------------
