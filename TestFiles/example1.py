import cantera
"""
This script demonstrates how to create a Cantera Solution object using the 'gri30.yaml' mechanism.

- `cantera.Solution`: Creates a Solution object representing a gas mixture defined by a mechanism file.
- To see all available attributes and methods for a Cantera Solution object, you can use Python's built-in `dir()` function:
    Example:
        print(dir(gas))
- For detailed documentation on each attribute and method, use the `help()` function:
    Example:
        help(gas)
- Commonly used properties include:
    - `species_names`: List of species in the mechanism.
    - `element_names`: List of elements in the mechanism.
    - `T`, `P`, `X`, `Y`: Temperature, pressure, mole fractions, and mass fractions.
    - `n_species`, `n_elements`: Number of species and elements.
    - `set_equivalence_ratio()`, `set_unnormalized_mass_fractions()`, etc.: Methods for setting composition.

Refer to the Cantera documentation for a comprehensive list of properties and methods:
https://cantera.org/documentation/
"""
gas = cantera.Solution("gri30.yaml")
#gas.TPX = 1200, 101325, 'CH4:1, O2:2' # Initial conditions, gets adabatic flame temperature of 3151,37k
gas.TPX = 298, 101325, 'CH4:1, O2:2'  # Initial conditions, gets adabatic flame temperature of 3052.05k
#gas()

gas.equilibrate('HP')
print(f"Adiabatic flame temperature: {gas.T:.2f} K")


