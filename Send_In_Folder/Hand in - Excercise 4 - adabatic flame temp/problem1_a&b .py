"""
This script was written for the purpose of assigments for the Combustion course on
Technical Univerity Munich.

Assigment: Excercise 4 - adabatic flame temp
Student: Herman Andersson

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import quad # Is this one used?

# Input fuel parameters
#x = int(input("Give the type of fuel C_xH_Y \n x = "))
#y = int(input("y = "))
#h_fuel = float(input("Give fuel heating value h_(f_i) [kJ/kmol]: "))

#----------------------------------
# Methane combustion parameters for testing code
x = 1
y = 4
h_fuel = -74831
#----------------------------------

# Enthalpy h_f,i values @ 298 K [kJ/kmol]
h_CO2 = -393546
h_H2O = -241845
h_N2 = 0
h_O2 = 0
h_CO = -110541
h_H2 = 0

# C_p,i values @ 1200 K [kJ/kmol-K]
C_p_CO2 = 56.21
C_p_H2O = 43.87
C_p_N2 = 33.71
C_p_CO = 34.148
C_p_H2 = 31.077

# Gas constant [J/mol-K]
R_u = 8.314
# Environment temperature [K]
T_0 = 298
# Starting temperature [K]
T_s = 1200

# C_p,i functions
def c_p_CO2(T):
    return R_u * (0.044536e2 + 0.03140168e-1*T - 0.12784105e-5*T**2 + 0.02393996e-8*T**3 - 0.16690333e-13*T**4)

def c_p_H2O(T):
    return R_u * (0.02672145e2 + 0.03056293e-1*T - 0.08730260e-5*T**2 + 0.12009964e-9*T**3 - 0.06391618e-13*T**4)

def c_p_N2(T):
    return R_u * (0.02926640e2 + 0.14879768e-2*T - 0.05684760e-5*T**2 + 0.10097038e-9*T**3 - 0.0675335e-13*T**4)

def c_p_O2(T):
    return R_u * (0.03697578e2 + 0.06135197e-2*T - 0.12588420e-6*T**2 + 0.01775281e-9*T**3 - 0.11364354e-14*T**4)

def c_p_CO(T):
    return R_u * (0.03025078e2 + 0.14426885e-2*T - 0.05630827e-5*T**2 + 0.10185813e-9*T**3 - 0.06910951e-13*T**4)

def c_p_H2(T):
    return R_u * (0.02991423e2 + 0.07000644e-2*T - 0.05633828e-6*T**2 - 0.09231578e-10*T**3 + 0.15827519e-14*T**4)

# Helper function to integrate c_p from T_0 to T # Is this one used?
def integrate_cp(cp_func, T):
    result, _ = quad(cp_func, T_0, T)
    return result

Phi = np.arange(0.1, 2.1, 0.1)
AFT = np.zeros(len(Phi))
AFT_O2 = np.zeros(len(Phi))

# ---------------------------------
# Problem 1.a // Calculate AFT with air
for index in range(len(Phi)):
    epsilon = 1
    T_s = 1200
    a = (x + y/4) / Phi[index]
    
    if Phi[index] <= 1:  # Fuel lean mixture
        # CxHy + a*(O_2 + 3.76*N_2) --> b*CO_2 + d*H_2O + f*O_2 + 3.76*a*N_2
        b = x
        d = y/2
        f = a - x - y/4
        H_R = h_fuel + a*h_O2 + 3.76*a*h_N2
        
        while epsilon > 0.01:
            def equation(T):
                H_P = (b*(h_CO2 + c_p_CO2(T_s)*(T - T_0)) + 
                       d*(h_H2O + c_p_H2O(T_s)*(T - T_0)) + 
                       f*(h_O2 + c_p_O2(T_s)*(T - T_0)) +
                       3.76*a*(h_N2 + c_p_N2(T_s)*(T - T_0)))
                return H_R - H_P

            
            T_new = fsolve(equation, T_s)[0]
            epsilon = abs(T_new - T_s)
            T_s = T_new
            
        AFT[index] = T_s
        
    else:  # Fuel rich mixture
        while epsilon > 0.01:
            # Water-gas shift equilibrium
            K_p = 10**(-2.4198 + 0.0003855*T_s + 2180.6/T_s)
            
            # Solve for b using equilibrium constant
            def solve_b(b_val):
                c_val = x - b_val
                d_val = 2*a - b_val - x
                e_val = -2*a + b_val + x + y/2
                if c_val <= 0 or d_val <= 0:
                    return float('inf')
                return K_p - (b_val * e_val) / (c_val * d_val)
            
            B = fsolve(solve_b, x/2)[0]
            C = x - B
            D = 2*a - B - x
            E = -2*a + B + x + y/2
            
            H_R = h_fuel + a*h_O2 + 3.76*a*h_N2
            
            def equation(T):
                H_P = (B*(h_CO2 + c_p_CO2(T_s)*(T - T_0)) + 
                       C*(h_CO + c_p_CO(T_s)*(T - T_0)) + 
                       D*(h_H2O + c_p_H2O(T_s)*(T - T_0)) + 
                       E*(h_H2 + c_p_H2(T_s)*(T - T_0)) +
                       3.76*a*(h_N2 + c_p_N2(T_s)*(T - T_0)))
                return H_R - H_P

            
            T_new = fsolve(equation, T_s)[0]
            epsilon = abs(T_new - T_s)
            T_s = T_new
            
        AFT[index] = T_s

# ---------------------------------
# Problem 1.b // Calculate AFT with pure O2 
for index in range(len(Phi)):
    epsilon = 1
    T_s = 1200
    a = (x + y/4) / Phi[index]
    
    if Phi[index] <= 1:  # Fuel lean mixture
        # CxHy + a*O_2 --> b*CO_2 + d*H_2O + f*O_2
        b = x
        d = y/2
        f = a - x - y/4
        H_R = h_fuel + a*h_O2
        
        while epsilon > 0.01:
            def equation(T):
                H_P = (b*(h_CO2 + c_p_CO2(T_s)*(T - T_0)) + 
                       d*(h_H2O + c_p_H2O(T_s)*(T - T_0)) + 
                       f*(h_O2 + c_p_O2(T_s)*(T - T_0)))
                return H_R - H_P

            
            T_new = fsolve(equation, T_s)[0]
            epsilon = abs(T_new - T_s)
            T_s = T_new
            
        AFT_O2[index] = T_s
        
    else:  # Fuel rich mixture
        while epsilon > 0.01:
            # Water-gas shift equilibrium
            K_p = 10**(-2.4198 + 0.0003855*T_s + 2180.6/T_s)
            
            # Solve for b using equilibrium constant
            def solve_b(b_val):
                c_val = x - b_val
                d_val = 2*a - b_val - x
                e_val = -2*a + b_val + x + y/2
                if c_val <= 0 or d_val <= 0:
                    return float('inf')
                return K_p - (b_val * e_val) / (c_val * d_val)
            
            B = fsolve(solve_b, x/2)[0]
            C = x - B
            D = 2*a - B - x
            E = -2*a + B + x + y/2
            
            H_R = h_fuel + a*h_O2
            
            def equation(T):
                H_P = (B*(h_CO2 + c_p_CO2(T_s)*(T - T_0)) + 
                       C*(h_CO + c_p_CO(T_s)*(T - T_0)) + 
                       D*(h_H2O + c_p_H2O(T_s)*(T - T_0)) + 
                       E*(h_H2 + c_p_H2(T_s)*(T - T_0)))
                return H_R - H_P

            
            T_new = fsolve(equation, T_s)[0]
            epsilon = abs(T_new - T_s)
            T_s = T_new
            
        AFT_O2[index] = T_s

# Increase the font size of the plot
plt.rcParams.update({
    'axes.titlesize': 24,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'axes.labelsize': 28 
})

# Plotting
plt.plot(Phi, AFT, 'b-', label='Air')
plt.axhline(y=max(AFT), color='r', linestyle='--', label='Air AFT Max')
plt.plot(Phi, AFT_O2, 'r-', label='Pure O₂')
plt.axhline(y=max(AFT_O2), color='b', linestyle='--', label='O2 AFT Max')
plt.title('Methane combustion', fontsize=40)
plt.xlabel('Equivalence ratio, φ', fontsize=30)
plt.ylabel('Adiabatic flame temperature [K]', fontsize=30)



plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

print("Calculation completed!")
#print(f"AFT with air: {AFT}")
#print(f"AFT with pure O2: {AFT_O2}")

# Should be 2200
print(f"Maximum AFT with pure O2: {np.max(AFT):.2f} K at phi = {Phi[np.argmax(AFT)]:.2f}")
# should be 4700
print(f"Maximum AFT with pure O2: {np.max(AFT_O2):.2f} K at phi = {Phi[np.argmax(AFT_O2)]:.2f}")

# Print PHI and AFT values on two lines
#print("PHI values:")
#print(Phi)
#print("AFT values:")
#print(AFT)