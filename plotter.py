'''
Header      : plotter.py
Created     : 06.08.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Plots two .csv datasets next to each other. Used to compare results with Hysys or Plus.
'''


import matplotlib.pyplot as plt
import pandas as pd


def plot_results(results: pd.DataFrame):
	"""
	Plots results of reactor simulation. Input DataFrame must be in format defined in ReactorModel.Simulation
	:param results: Tabular results of PFReactor.Simulation
	"""
	fig, axes = plt.subplots(2, 1)
	axes[0].set_title('Composition and Temperature of Reaction Mixture vs. Time')
	axes[0].grid()
	axes[0].set_ylabel('Molar Fraction, [mol. frac.]')
	for param in results.columns[0:-1]:
		if 'rate' not in param and 'FLMASS' not in param and 'FLMOL' not in param and param != 't [s]':
			axes[0].plot(results.index, results[param], label= param)
	axes[0].legend()
	axes[1].set_ylabel('Temperature, [K]')
	axes[1].set_xlabel('Length, [m]')
	axes[1].plot(results.index, results['T [K]'])
	plt.grid()
	plt.show()
	return 1


### WARNING! Needs to be wrapped into function
def printplot():
	main_plot_filename = 'log.csv'
	secondary_plot_filename = 'log_hysys_dante_1000_steps.csv' #'log_solver_runge-kutta-1e-2.csv'

	main_plot_data = pd.read_csv(main_plot_filename)
	# secondary_plot_data = pd.read_csv(secondary_plot_filename)
	secondary_plot_data = pd.read_csv(secondary_plot_filename, sep= ';')

	fig, axes = plt.subplots(3, 1)
	axes[0].set_title('Composition, Temperature of Reaction Mixture and Reaction Rates vs. Length')
	axes[0].grid()
	axes[0].set_ylabel('Molar Fraction, [mol. frac.]')
	for param in main_plot_data.columns[1:-1]:
		if 'rate' not in param and 'FLMASS' not in param and 'FLMOL' not in param and param != 't [s]':
			axes[0].plot(main_plot_data['l [m]'], main_plot_data[param], label=param + '-solver', color= 'red', linewidth= 0.5)
			axes[0].plot(secondary_plot_data['l [m]'], secondary_plot_data[param], label=param + '-hysys',
						 color= 'blue', linewidth= 0.5, linestyle= '--')

	# axes[0].legend()
	axes[1].set_ylabel('Temperature, [K]')
	axes[1].set_xlabel('Length, [m]')
	axes[1].grid()
	axes[1].plot(main_plot_data['l [m]'], main_plot_data['T [K]'], color= 'red', linewidth= 0.5)
	axes[1].plot(secondary_plot_data['l [m]'], secondary_plot_data['T [K]'], color= 'blue', linewidth= 0.5, linestyle= '--')
	axes[2].set_ylabel('Reaction Rates, [kgmol/(m3*s)]')
	for param in main_plot_data.columns[0:-1]:
		if 'rate' in param:
			axes[2].plot(main_plot_data['l [m]'], main_plot_data[param], label=param + '-solver', color= 'red', linewidth= 0.5)
			axes[2].plot(secondary_plot_data['l [m]'], secondary_plot_data[param], label=param + '-hysys',
						 color= 'blue', linewidth= 0.5, linestyle= '--')
	axes[2].legend()

	plt.grid()
	plt.show()
	return 1

# printplot()