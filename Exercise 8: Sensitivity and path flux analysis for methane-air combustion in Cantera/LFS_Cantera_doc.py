import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['figure.constrained_layout.use'] = True

# Define the gas mixture and kinetics used to compute mixture properties
# In this case, we are using the GRI 3.0 model with methane as the fuel
mech_file = "gri30.yaml"
fuel_comp = {"CH4": 1.0}

# Inlet temperature in kelvin and inlet pressure in pascals
# In this case we are setting the inlet T and P to room temperature conditions
To = 300
Po = ct.one_atm

# Domain width in metres
width = 0.014

# Create the object representing the gas and set its state to match the inlet conditions
gas = ct.Solution(mech_file)
gas.TP = To, Po
gas.set_equivalence_ratio(1.0, fuel_comp, {"O2": 1.0, "N2": 3.76})

flame = ct.FreeFlame(gas, width=width)
flame.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
flame.solve(loglevel=1, auto=True)

# Create a DataFrame to store sensitivity-analysis data
sens = pd.DataFrame(index=gas.reaction_equations(), columns=["sensitivity"])

# Use the adjoint method to calculate sensitivities
sens.sensitivity = flame.get_flame_speed_reaction_sensitivities()

# Show the first 10 sensitivities
sens.head(10)

# Sort the sensitivities in order of descending magnitude
sens = sens.iloc[(-sens['sensitivity'].abs()).argsort()]

fig, ax = plt.subplots()

# Reaction mechanisms can contains thousands of elementary steps. Limit the plot
# to the top 15
sens.head(15).plot.barh(ax=ax, legend=None)

ax.invert_yaxis()  # put the largest sensitivity on top
ax.set_title("Sensitivities for GRI 3.0")
ax.set_xlabel(r"Sensitivity: $\frac{\partial\:\ln S_u}{\partial\:\ln k}$")
ax.grid(axis='x')
plt.show()