import pandas as pd
import numpy as np

# needed to save the figures without displaying them
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting

import matplotlib.pyplot as plt
import cantera as ct

# -------------------- Experimental Data import Section --------------------------------
# Import experimental flame speed data
exp_const_temp = pd.read_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/FlameSpeedsExperimental_constant temperature_298K.csv')
exp_const_P = pd.read_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/FlameSpeedsExperimental_constant pressure_1atm.csv')

# Only take coloumns with relevant data
exp_const_temp = exp_const_temp.iloc[:, [0, 1, 2, 3, 8, 9, 12, 13]]
exp_const_P = exp_const_P.iloc[:, :-3]


# -------------------- Cantera Calculation Section --------------------------------
# Simulation Parameters
pressure = [101325, 202650, 506635, 2026500]  # pressure in Pa]
Tin_const = 298.0  # unburned gas temperature in K
Tin_var = [358.0, 393.0, 428.0]  # unburned gas temperature in K
width = 0.03  # width of the flame in m

# Solution object used t compute mixture properties, 
# set to to the state of upstream fuel-air mixture
gas = ct.Solution("gri30.yaml")  

phi = np.linspace(0.5, 1.5, 21)  # equivalence ratio range
loglevel = 1  # amount of diagnostic outptut (0 to 8)
fin_P_Result=np.zeros((len(pressure), len(phi))) # Empty array to store results
fin_Temp_Result=np.zeros((len(Tin_var), len(phi))) # Empty array to store results


def calculate_flame_speed(param_list, phi, gas, loglevel, fixed_value, is_pressure_variable, result_array, width):
    for i in range(len(param_list)):
        for j in range(len(phi)):
            print("-"*50 + "*"*50 + "-"*50)
            print("-"*50 + "*"*50 + "-"*50)
            print("-"*50 + "Current iteration: " + str(i) + ", " + str(j) + "-"*50)
            print("-"*50 + "*"*50 + "-"*50)
            print("-"*50 + "*"*50 + "-"*50)
            # Set equivalence ratio
            gas.set_equivalence_ratio(phi[j], 'CH4', 'O2:2, N2:7.52')
            
            # Assign T and P based on which is being varied
            if is_pressure_variable:
                gas.TP = fixed_value, param_list[i]  # Temp is fixed, pressure varies
            else:
                gas.TP = param_list[i], fixed_value  # Pressure is fixed, temp varies

            # Set up flame object
            flame = ct.FreeFlame(gas, width=width)
            flame.set_refine_criteria(ratio=2.0, slope=0.06, curve=0.12)
            flame.solve(loglevel=loglevel, auto=True)

            # Store flame speed result in cm/s
            result_array[i][j] = flame.velocity[0] * 100


# --- Case 1: Varying Pressure ---
calculate_flame_speed(pressure, phi, gas, loglevel, Tin_const, True, fin_P_Result, width)
# --- Case 2: Varying Temperature ---
calculate_flame_speed(Tin_var, phi, gas, loglevel, pressure[0], False, fin_Temp_Result, width)


# -------------------- Plot Section --------------------------------

# Marker styles and colors to differentiate datasets
markers = ['o', 's', '^', 'D', 'v', '*', 'x', 'P']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:gray']


fig, axs = plt.subplots(2, 1, figsize=(10, 12))

def data_to_plot(exp_data, markers, colors, ax, labels):
    num_columns = exp_data.shape[1]
    for i in range(0, num_columns, 2):
        x = exp_data.iloc[:, i]
        y = exp_data.iloc[:, i + 1]
        idx = i // 2
        marker = markers[idx % len(markers)]
        color = colors[idx % len(colors)]
        label = labels[idx % len(labels)] + " (Exp)" 
        ax.plot(x, y, linestyle='None', marker=marker, color=color, label=label)


# Plot first dataset (top subplot)
labels_temp = ["358 k", "393 k", "428 k"]
data_to_plot(exp_const_P, markers, colors, axs[0], labels_temp)
for i in range(len(Tin_var)):
    plots = axs[0].plot(phi, fin_Temp_Result[i], marker=markers[i], linestyle='-', color=colors[i], label= labels_temp[i] + ' (Cantera)')
axs[0].set_title('Flame Speed at Constant Pressure (1 atm)')
axs[0].set_xlabel('Equivalence Ratio ϕ')
axs[0].set_ylabel('Laminar Flame Speed [cm/s]')
axs[0].legend()
axs[0].grid(True)

# Plot second dataset (bottom subplot)
labels_presure = ["1 atm", "2 atm", "5 atm", "20 atm"]
data_to_plot(exp_const_temp, markers, colors, axs[1], labels_presure)
for i in range(len(pressure)):
    plots = axs[1].plot(phi, fin_P_Result[i], marker=markers[i], linestyle='-', color=colors[i], label= labels_presure[i] + ' (Cantera)')
axs[1].set_title('Flame Speed at Constant Temperature (298 K)')
axs[1].set_xlabel('Equivalence Ratio ϕ')
axs[1].set_ylabel('Laminar Flame Speed [cm/s]')
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.show()

# -------------------- Save Result Arrays to CSV with Labels --------------------

# Save pressure-varying results with pressure as first column and phi as column headers
df_P = pd.DataFrame(fin_P_Result, columns=[f'phi={round(p, 3)}' for p in phi])
df_P.insert(0, "Pressure_Pa", pressure)
df_P.to_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/fin_P_Result.csv', index=False)

# Save temperature-varying results with temperature as first column and phi as column headers
df_T = pd.DataFrame(fin_Temp_Result, columns=[f'phi={round(p, 3)}' for p in phi])
df_T.insert(0, "Temp_K", Tin_var)
df_T.to_csv('/home/brukare/Combustion/Excercise 7 - laminar flame speed/fin_Temp_Result.csv', index=False)


# Save the figure, won't display it, need to activate Agg backend under imports
plt.savefig('/home/brukare/Combustion/Excercise 7 - laminar flame speed/FlameSpeed_Results.png', dpi=300)
