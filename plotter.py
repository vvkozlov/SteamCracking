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


def plotlog(log_name: str):
	"""
	Plots data from log file in .csv format created in chemistry.simulate() method.
	One should pay close attention to formatting of log file because column names are used to differentiate data.
	Note that log file has to be in same directory as plotter.py script

	:param log_name: [str] Name of log file (with extension) to be printed
	:return: nothing
	"""

	"""Read log file to pd.DataFrame"""
	try:
		log_data = pd.read_csv(log_name, sep= ',')  # Try with first column separator
		separator_check = log_data['l [m]']  # Check validity of .csv columns separator
	except KeyError:
		log_data = pd.read_csv(log_name, sep=';')  # Try with second column separator
		separator_check = log_data['l [m]']  # Check validity of .csv columns separator
	except:
		print(f'ERROR! Unable to locate logfile named {log_name}. Check filename spelling')
		return 1
	"""Extract data headers"""
	log_data_headers = list(log_data.columns)[1:]  # List of all headers to choose from
	"""Set names for process parameters columns"""
	temperature_header = 'T [K]'
	molflow_header = 'FLMOL [kgmol/hr]'
	residence_time_header = 't [s]'
	process_data_headers = list([temperature_header, molflow_header, residence_time_header])
	"""Prepare data headers to differentiate data"""
	rxn_data_headers = [x for x in log_data_headers if '-->' in x]  # Stream composition data
	comp_data_headers = [x for x in log_data_headers if x not in process_data_headers
						 and x not in rxn_data_headers]  # Reaction rates data
	length_steps = log_data['l [m]']  # List of integration lengths to use as indexes
	"""Layout of graphs"""
	comp_plot = plt.subplot2grid((8, 8), (0, 0), rowspan= 4, colspan= 4)
	plt.subplots_adjust(wspace= 0.7)
	rxn_plot = plt.subplot2grid((8, 8), (0, 4), rowspan= 4, colspan= 4)
	plt.subplots_adjust(hspace= 1.5)
	temperature_plot = plt.subplot2grid((8, 8), (4, 0), rowspan= 2, colspan= 8)
	molflow_plot = plt.subplot2grid((8, 8), (6, 0), rowspan= 2, colspan= 8)
	"""Set axes limits"""
	# At the moment auto limits looks just fine
	"""Set graphs labels"""
	comp_plot.set_ylabel('Composition [mol. fract.]')
	comp_plot.set_xlabel('Length [m]')
	rxn_plot.set_ylabel('Reaction Rate [kgmol/(m3*s)]')
	rxn_plot.set_xlabel('Length [m]')
	temperature_plot.set_ylabel('Temperature [K]')
	molflow_plot.set_ylabel('Molar Flow [kgmol/hr]')
	molflow_plot.set_xlabel('Length [m]')
	"""Print data"""
	comp_plot.plot(length_steps, log_data[comp_data_headers])
	rxn_plot.plot(length_steps, log_data[rxn_data_headers])
	temperature_plot.plot(length_steps, log_data['T [K]'])
	try:
		molflow_plot.plot(length_steps, log_data['FLMOL [kgmol/hr]'])
	except:
		molflow_plot.text(0.1, 0.1, 'Molar Flow Data is not availible')
	plt.show()

plotlog('log.csv')

### WARNING! Does not work properly
def compare_logs():
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