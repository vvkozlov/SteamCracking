import matplotlib.pyplot as plt
import os
import pandas as pd

main_plot_filename = 'log_solver_chlorine.csv'
secondary_plot_filename = 'log_hysys_chlorine.csv'

main_plot_data = pd.read_csv(main_plot_filename)
secondary_plot_data = pd.read_csv(secondary_plot_filename, sep= ';')

fig, axes = plt.subplots(3, 1)
axes[0].set_title('Composition, Temperature of Reaction Mixture and Reaction Rates vs. Length')
axes[0].grid()
axes[0].set_ylabel('Molar Fraction, [mol. frac.]')
for param in main_plot_data.columns[1:-1]:
	if 'rate' not in param and 'FLMASS' not in param and 'FLMOL' not in param:
		axes[0].plot(main_plot_data['l'], main_plot_data[param], label=param, color= 'red', linewidth= 0.5)
		axes[0].plot(secondary_plot_data['l'], secondary_plot_data[param], label=param,
					 color= 'blue', linewidth= 0.5, linestyle= '--')

axes[0].legend()
axes[1].set_ylabel('Temperature, [K]')
axes[1].set_xlabel('Length, [m]')
axes[1].grid()
axes[1].plot(main_plot_data['l'], main_plot_data['T'], color= 'red', linewidth= 0.5)
axes[1].plot(secondary_plot_data['l'], secondary_plot_data['T'], color= 'blue', linewidth= 0.5, linestyle= '--')
axes[2].set_ylabel('Reaction Rates, [kgmol/(m3*s)]')
for param in main_plot_data.columns[0:-1]:
	if 'rate' in param:
		axes[2].plot(main_plot_data['l'], main_plot_data[param], label=param, color= 'red', linewidth= 0.5)
		axes[2].plot(secondary_plot_data['l'], secondary_plot_data[param], label=param,
					 color= 'blue', linewidth= 0.5, linestyle= '--')
axes[2].set_ylim([0, 0.05])

plt.grid()
plt.show()

