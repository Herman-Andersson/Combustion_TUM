import numpy as np
import matplotlib.pyplot as plt

def calculate_cp(T, coeffs):
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
    # For methane: x=1, y=4
    x, y = 1, 4
    a = (x + y/4) / phi  # Stoichiometric oxygen requirement
    
    # Water-gas equilibrium constant
    Kp = water_gas_equilibrium_kp(T)
    
    if phi <= 1.0:  # Lean mixture
        # CH4 + a*O2 + 3.76*a*N2 -> b*CO2 + d*H2O + f*O2 + 3.76*a*N2
        b = x  # All carbon goes to CO2
        d = y/2  # All hydrogen goes to H2O
        f = a - b - d/2  # Remaining oxygen
        c = 0  # No CO
        e = 0  # No H2
    else:  # Rich mixture
        # Need to solve for b using water-gas equilibrium
        # From mass balance equations:
        # c = x - b  (equation 1)
        # d = 2*a - b - x  (equation 2)
        # e = -2*a + b + x + y/2  (equation 3)
        # Kp = (b*e)/(c*d)  (equation 4)
        
        # Substituting equations 1, 2, 3 into 4:
        # Kp = [b * (-2*a + b + x + y/2)] / [(x - b) * (2*a - b - x)]
        
        # This gives a quadratic equation in b:
        # Kp * (x - b) * (2*a - b - x) = b * (-2*a + b + x + y/2)
        # Expanding and rearranging: A*b^2 + B*b + C = 0
        
        A = Kp + 1
        B = -Kp*(3*x + 2*a) - (-2*a + x + y/2)
        C = Kp * x * (2*a - x)
        
        # Solve quadratic equation
        discriminant = B**2 - 4*A*C
        if discriminant < 0:
            # Fallback to complete combustion
            b = min(x, 2*a)
        else:
            b1 = (-B + np.sqrt(discriminant)) / (2*A)
            b2 = (-B - np.sqrt(discriminant)) / (2*A)
            
            # Choose physically meaningful root (0 <= b <= x)
            if 0 <= b1 <= x:
                b = b1
            elif 0 <= b2 <= x:
                b = b2
            else:
                b = min(x, 2*a)
        
        # Calculate other species
        c = x - b
        d = 2*a - b - x
        e = -2*a + b + x + y/2
        f = 0  # No excess oxygen in rich mixture
        
        # Ensure non-negative values
        c = max(0, c)
        d = max(0, d)
        e = max(0, e)
    
    return b, c, d, e, f, 3.76*a

def calculate_aft(phi_values):
    """Calculate adiabatic flame temperature for given equivalence ratios"""
    
    # Thermodynamic data
    h_f = {  # Formation enthalpies at 298K (kJ/kmol)
        'CH4': -74831,
        'CO2': -393546,
        'H2O': -241845,
        'N2': 0,
        'O2': 0,
        'CO': -110541,
        'H2': 0
    }
    
    # Polynomial coefficients for cp calculation (1000-5000K)
    coeffs = {
        'CO2': [0.044536E+02, 0.03140168E-01, -0.12784105E-05, 0.02393996E-08, -0.16690333E-13],
        'H2O': [0.02672145E+02, 0.03056293E-01, -0.08730260E-05, 0.12009964E-09, -0.06391618E-13],
        'N2':  [0.02926640E+02, 0.14879768E-02, -0.05684760E-05, 0.10097038E-09, -0.0675335E-13],
        'O2':  [0.03697578E+02, 0.06135197E-02, -0.12588420E-06, 0.01775281E-09, -0.11364354E-14],
        'CO':  [0.03025078E+02, 0.14426885E-02, -0.05630827E-05, 0.10185813E-09, -0.06910951E-13],
        'H2':  [0.02991423E+02, 0.07000644E-02, -0.05633828E-06, -0.09231578E-10, 0.15827519E-14]
    }
    
    T_initial = 298.15  # K
    results = []
    
    for phi in phi_values:
        print(f"Calculating AFT for φ = {phi:.2f}")
        
        # Initial guess for AFT
        T_aft = 1200.0  # K
        tolerance = 1.0  # K
        max_iterations = 100
        
        for iteration in range(max_iterations):
            T_old = T_aft
            
            # Solve species balance
            n_CO2, n_CO, n_H2O, n_H2, n_O2, n_N2 = solve_species_balance(phi, T_aft)
            
            # Calculate enthalpy of reactants (at 298K)
            # CH4 + a*O2 + 3.76*a*N2
            x, y = 1, 4
            a = (x + y/4) / phi
            
            H_reactants = (1 * h_f['CH4'] + a * h_f['O2'] + 3.76*a * h_f['N2'])
            
            # Calculate enthalpy of products at AFT
            H_products = 0
            species_data = [
                ('CO2', n_CO2), ('CO', n_CO), ('H2O', n_H2O), 
                ('H2', n_H2), ('O2', n_O2), ('N2', n_N2)
            ]
            
            for species, n_moles in species_data:
                if n_moles > 0:
                    # Formation enthalpy
                    H_products += n_moles * h_f[species]
                    
                    # Sensible enthalpy (integrate cp from 298K to T_aft)
                    if species in coeffs:
                        cp_coeffs = coeffs[species]
                        # Approximate integration using average cp
                        T_avg = (T_initial + T_aft) / 2
                        cp_avg = calculate_cp(T_avg, cp_coeffs)
                        H_products += n_moles * cp_avg * (T_aft - T_initial)
            
            # Energy balance: H_reactants = H_products
            # Adjust T_aft to satisfy energy balance
            if iteration > 0:
                # Calculate total heat capacity of products
                cp_total = 0
                for species, n_moles in species_data:
                    if n_moles > 0 and species in coeffs:
                        cp_total += n_moles * calculate_cp(T_aft, coeffs[species])
                
                if cp_total > 0:
                    # Newton-Raphson type adjustment
                    dH = H_reactants - H_products
                    T_aft = T_aft + dH / cp_total
                else:
                    T_aft = T_aft + 50  # Small increment if cp_total is zero
            
            # Check convergence
            if abs(T_aft - T_old) < tolerance:
                print(f"  Converged after {iteration+1} iterations: T_AFT = {T_aft:.2f} K")
                break
            
            # Prevent unrealistic temperatures
            T_aft = max(300, min(4000, T_aft))
        
        else:
            print(f"  Did not converge after {max_iterations} iterations: T_AFT = {T_aft:.2f} K")
        
        results.append({'phi': phi, 'T_AFT': T_aft})
    
    return results

# Main execution
if __name__ == "__main__":
    # Define equivalence ratio range
    phi_values = np.arange(0.1, 2.1, 0.1)
    
    print("Calculating Adiabatic Flame Temperature for Methane-Air Combustion")
    print("=" * 60)
    
    # Calculate AFT for each equivalence ratio
    results = calculate_aft(phi_values)
    
    # Extract phi and T_AFT values for plotting
    phi_array = np.array([result['phi'] for result in results])
    T_AFT_array = np.array([result['T_AFT'] for result in results])
    
    # Display results
    print("\nResults:")
    print("-" * 30)
    for result in results:
        phi = result['phi']
        T_AFT = result['T_AFT']
        print(f"φ = {phi:.2f}, T_AFT = {T_AFT:.2f} K ({T_AFT-273.15:.2f} °C)")
    
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(phi_array, T_AFT_array, 'b-o', linewidth=2, markersize=6)
    plt.xlabel('Equivalence Ratio (φ)')
    plt.ylabel('Adiabatic Flame Temperature (K)')
    plt.title('Adiabatic Flame Temperature vs Equivalence Ratio\nMethane-Air Combustion')
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 2.1)
    plt.tight_layout()
    plt.show()
    
    # Find maximum AFT
    max_aft_idx = np.argmax(T_AFT_array)
    max_aft_phi = phi_array[max_aft_idx]
    max_aft_temp = T_AFT_array[max_aft_idx]
    
    print(f"\nMaximum AFT: {max_aft_temp:.2f} K ({max_aft_temp-273.15:.2f} °C) at φ = {max_aft_phi:.2f}")