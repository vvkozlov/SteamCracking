'''
Header      : plotter.py
Created     : 06.08.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Plots two .csv datasets next to each other. Used to compare results with Hysys or Plus.
'''

import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import scienceplots

# plt.rcParams.update({'figure.dpi': '100'})


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
    # log_path = os.path.join(os.getcwd(), 'log\{}.csv'.format(log_name))
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
                         and x not in rxn_data_headers
                         and x not in ['Length step', 'ODE convergence steps', 'MBAL convergence steps']]  # Reaction rates data
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
    """Add legend"""
    rxn_plot.legend(labels= rxn_data_headers, title= 'huh?', loc='upper left', bbox_to_anchor=(1, 0.5))
    plt.show()

def plothist(history: pd.DataFrame):
    """
    Plots data from cal_hist in pd.DataFrame format created in chemistry.simulate() method.

    :param history: [pd.DataFrame] Calculations history variable to be printed
    :return: nothing
    """
    """Extract data headers"""
    log_data_headers = list(history.columns)[1:]  # List of all headers to choose from
    """Set names for process parameters columns"""
    temperature_header = 'T [K]'
    molflow_header = 'FLMOL [kgmol/hr]'
    residence_time_header = 't [s]'
    process_data_headers = list([temperature_header, molflow_header, residence_time_header])
    """Prepare data headers to differentiate data"""
    rxn_data_headers = [x for x in log_data_headers if '-->' in x]  # Stream composition data
    comp_data_headers = [x for x in log_data_headers if x not in process_data_headers
                         and x not in rxn_data_headers
                         and x not in ['Length step', 'ODE convergence steps', 'MBAL convergence steps']]  # Reaction rates data
    length_steps = history['Length step']  # List of integration lengths to use as indexes
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
    comp_plot.plot(length_steps, history[comp_data_headers])
    rxn_plot.plot(length_steps, history[rxn_data_headers])
    temperature_plot.plot(length_steps, history['T [K]'])
    try:
        molflow_plot.plot(length_steps, history['FLMOL [kgmol/hr]'])
    except:
        molflow_plot.text(0.1, 0.1, 'Molar Flow Data is not availible')
    """Add legend"""
    rxn_plot.legend(labels= rxn_data_headers, title= 'huh?', loc='upper left', bbox_to_anchor=(1, 0.5))
    plt.show()


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

class Logger(object):
    """
    Lumberjack class - duplicates sys.stdout to a log file and it okay
    source: /questions/49207/how-do-i-duplicate-sysstdout-to-a-log-file-in-python/348798#348798
    """
    def __init__(self, filename="console.dat", mode="w", buff=1):
        self.stdout = sys.stdout
        self.file = open(filename, mode, buff)
        sys.stdout = self

    def __del__(self):
        self.close()

    def __enter__(self):
        pass

    def __exit__(self, *args):
        pass

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()
        os.fsync(self.file.fileno())

    def close(self):
        if self.stdout != None:
            sys.stdout = self.stdout
            self.stdout = None
        if self.file != None:
            self.file.close()
            self.file = None


def plot_profiles(log_name: str):
    """Read log file to pd.DataFrame"""
    # log_path = os.path.join(os.getcwd(), 'log\{}.csv'.format(log_name))
    try:
        log_data = pd.read_csv(log_name, sep=',')  # Try with first column separator
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
    comp_data_headers = [x for x in log_data_headers if x not in process_data_headers and x != 'LTE'
                         and x not in rxn_data_headers]  # Reaction rates data
    #reduce headers
    comp_data_headers = comp_data_headers[:6]
    rxn_data_headers = rxn_data_headers[:6]
    length_steps = log_data['l [m]']  # List of integration lengths to use as indexes
    """Layout of graphs"""
    plt.style.use(['science', 'no-latex', 'grid', 'high-vis'])
    # f1 = plt.figure(figsize=(10, 4))
    f1 = plt.figure(figsize=(4, 4))
    comp_plot = plt.subplot2grid((10, 10), (0, 0), rowspan= 10, colspan= 10)
    """Set axes limits"""
    # At the moment auto limits looks just fine
    """Set graphs labels"""
    comp_plot.set_ylabel('Composition [mol. fract.]')
    comp_plot.set_xlabel('Length [m]')
    """Print data"""
    comp_plot.plot(length_steps, log_data[comp_data_headers])
    # comp_plot.legend(labels= comp_data_headers, title= 'Components', loc='center right', borderaxespad=-10)  #, ncol=2)
    comp_plot.legend(labels= comp_data_headers, title= 'Components', loc='upper right', borderaxespad=1, prop={'size': 6})  #, ncol=2)
    # plt.show()
    """Plotting rxn rates"""
    f2 = plt.figure(figsize=(4, 4))
    rxn_plot = plt.subplot2grid((8, 8), (0, 0), rowspan= 8, colspan= 8)
    rxn_plot.set_ylabel('Reaction Rate [kgmol/(m3*s)]')
    rxn_plot.set_xlabel('Length [m]')
    rxn_plot.plot(length_steps, log_data[rxn_data_headers])
    rxn_plot.legend(labels= rxn_data_headers, title= 'Reaction Rates', loc='upper right', borderaxespad=1, prop={'size': 6})
    f3 = plt.figure(figsize=(10, 4))
    temperature_plot = plt.subplot2grid((8, 8), (0, 0), rowspan= 8, colspan= 8)
    # temperature_plot = rxn_plot.twinx()
    temperature_plot.set_ylabel('Temperature [K]')
    temperature_plot.set_xlabel('Length [m]')
    temperature_plot.plot(length_steps, log_data[temperature_header])
    temperature_plot.legend(labels= temperature_header, title= 'Temperature', loc='upper right')
    # temperature_plot.set_ylim(bottom=0, top= 1400)
    # plt.ylabel('khui')
    plt.show()


def plot_compare(log_name_1: str, log_name_2: str, header: str, labels: list[str], marker: str):
    """Import first log
    :param marker allows to select dots or line for data
    """
    try:
        log_data_1 = pd.read_csv(log_name_1, sep=',')  # Try with first column separator
        separator_check = log_data_1['l [m]']  # Check validity of .csv columns separator
    except KeyError:
        log_data_1 = pd.read_csv(log_name_1, sep=';')  # Try with second column separator
        separator_check = log_data_1['l [m]']  # Check validity of .csv columns separator
    except:
        print(f'ERROR! Unable to locate logfile named {log_name_1}. Check filename spelling')
        return 1
    """Import second log"""
    try:
        log_data_2 = pd.read_csv(log_name_2, sep=',')  # Try with first column separator
        separator_check = log_data_2['l [m]']  # Check validity of .csv columns separator
    except KeyError:
        log_data_2 = pd.read_csv(log_name_2, sep=';')  # Try with second column separator
        separator_check = log_data_2['l [m]']  # Check validity of .csv columns separator
    except:
        print(f'ERROR! Unable to locate logfile named {log_name_2}. Check filename spelling')
        return 1
    """Extract data headers"""
    log_data_headers = list(log_data_1.columns)[1:]  # List of all headers to choose from
    temperature_header = 'T [K]'
    molflow_header = 'FLMOL [kgmol/hr]'
    residence_time_header = 't [s]'
    process_data_headers = list([temperature_header, molflow_header, residence_time_header])
    """Prepare data headers to differentiate data"""
    rxn_data_headers = [x for x in log_data_headers if '-->' in x]  # Stream composition data
    comp_data_headers = [x for x in log_data_headers if x not in process_data_headers and x != 'LTE'
                         and x not in rxn_data_headers]  # Reaction rates data
    """Export steps"""
    length_steps_1 = log_data_1['l [m]']  # List of integration lengths to use as indexes
    length_steps_2 = log_data_2['l [m]']  # List of integration lengths to use as indexes
    """Layout of graphs"""
    plt.style.use(['science', 'no-latex', 'grid', 'high-vis'])
    # f1 = plt.figure(figsize=(10, 4))
    f1 = plt.figure(figsize=(4, 4))
    # plot_to_compare = plt.subplot2grid((8, 8), (0, 0), rowspan=8, colspan=8)
    plot_to_compare = plt.subplot2grid((10, 12), (0, 1), rowspan=10, colspan=11)
    """Set axes limits"""
    # At the moment auto limits looks just fine
    """Set graphs labels"""
    custom_ylabel = ''
    if header in process_data_headers:
        if header == temperature_header:
            custom_ylabel = 'Temperature [K]'
        else:
            pass
    elif header[0] in comp_data_headers:  # works only for single header rn
        custom_ylabel = '{} Content [mol. fract.]'.format(header[0])
    elif header in rxn_data_headers:
        custom_ylabel = '{} Reaction Rate [kgmol/(m3*s)]'.format(header)
    elif header == 'LTE':
        custom_ylabel = 'Local Truncation Error Estimation'
        lte_median = (max(log_data_1[header]) + min(log_data_1[header])) / 2
        log_length_1 = len(log_data_1[header])
        points_above_median_1 = sum(1 for item in log_data_1[header] if item >= lte_median)
        ratio_above_median_1 = points_above_median_1 / log_length_1
        log_length_2 = len(log_data_2[header])
        print('Total points for log 1\t-\t', log_length_1)
        print('Median LTE for log 1\t-\t', lte_median)
        print('Points with LTE above median for log 1\t-\t', points_above_median_1)
        print('Ratio of points with LTE above median for log 1\t-\t', ratio_above_median_1)
        print('Total points for log 2\t-\t', log_length_2)
    elif header == 'step':
        custom_ylabel = 'Integration Step Size [m]'
        step_size_1 = list([])
        for i in range(len(length_steps_1) - 1):
            step_size_1.append(length_steps_1[i+1] - length_steps_1[i])
        step_size_1.append(step_size_1[-1])
        log_data_1['step'] = step_size_1
        step_size_2 = list([])
        for i in range(len(length_steps_2) - 1):
            step_size_2.append(length_steps_2[i + 1] - length_steps_2[i])
        step_size_2.append(step_size_2[-1])
        log_data_2['step'] = step_size_2
    elif header == 'stepnum':
        custom_ylabel = 'Number of Integration Steps'
        num_steps_1 = list([])
        steps_counter_1 = 0
        for i in range(len(length_steps_1)):
            steps_counter_1 += 1
            num_steps_1.append(steps_counter_1)
        log_data_1['stepnum'] = num_steps_1
        num_steps_2 = list([])
        steps_counter_2 = 0
        for i in range(len(length_steps_2)):
            steps_counter_2 += 1
            num_steps_2.append(steps_counter_2)
        log_data_2['stepnum'] = num_steps_2
    else:
        pass
    plot_to_compare.set_ylabel(custom_ylabel)
    plot_to_compare.set_xlabel('Length [m]')
    """Print data"""
    if marker == 'points':
        plot_to_compare.plot(length_steps_1, log_data_1[header], 'or', length_steps_2, log_data_2[header], 'ob', ms=0.2)
    elif marker == 'line':
        plot_to_compare.plot(length_steps_1, log_data_1[header], length_steps_2, log_data_2[header])
    else:
        print('ERROR! Select marker option')
        pass
    plot_to_compare.legend(labels=[labels[0], labels[1]], title='Cases', loc='lower right', borderaxespad=1, prop={'size': 6})
    plt.show()


def plot_selected(log_name: str, headers: list[str]):
    """Read log file to pd.DataFrame"""
    # log_path = os.path.join(os.getcwd(), 'log\{}.csv'.format(log_name))
    try:
        log_data = pd.read_csv(log_name, sep=',')  # Try with first column separator
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
    comp_data_headers = [x for x in log_data_headers if x not in process_data_headers and x != 'LTE'
                         and x not in rxn_data_headers]  # Reaction rates data
    #reduce headers
    length_steps = log_data['l [m]']  # List of integration lengths to use as indexes
    """Layout of graphs"""
    plt.style.use(['science', 'no-latex', 'grid', 'high-vis'])
    # f1 = plt.figure(figsize=(10, 4))
    f1 = plt.figure(figsize=(4, 4))
    # plot_selected = plt.subplot2grid((8, 8), (0, 0), rowspan=8, colspan=8)
    plot_selected = plt.subplot2grid((8, 12), (0, 1), rowspan=8, colspan=11)
    """Set axes limits"""
    # At the moment auto limits looks just fine
    """Set graphs labels"""
    custom_ylabel = ''
    if headers in process_data_headers:
        if len(headers) == 1 and headers[0] == temperature_header:
            custom_ylabel = 'Temperature [K]'
        else:
            pass
    elif headers[0] in comp_data_headers: # cheating here
        custom_ylabel = 'Composition [mol. fract.]'
    elif headers[0] in rxn_data_headers:
        custom_ylabel = 'Reaction Rates [kgmol/(m3*s)]'
    elif len(headers) == 1 and headers[0] == 'LTE':
        custom_ylabel = 'Local Truncation Error Estimation'
        lte_median = (max(log_data[headers[0]]) + min(log_data[headers[0]])) / 2
        log_length = len(log_data[headers[0]])
        points_above_median = sum(1 for item in log_data[headers[0]] if item >= lte_median)
        ratio_above_median = points_above_median/ log_length
        print('Total points for log\t-\t', log_length)
        print('Median LTE for log\t-\t', lte_median)
        print('Points with LTE above median for log\t-\t', points_above_median)
        print('Ratio of points with LTE above median for log\t-\t', ratio_above_median)
    elif len(headers) == 1 and headers[0] == 'step':
        custom_ylabel = 'Integration Step Size [m]'
        step_size = list([])
        for i in range(len(length_steps) - 1):
            step_size.append(length_steps[i + 1] - length_steps[i])
        step_size.append(step_size[-1])
        log_data['step'] = step_size
        step_size_2 = list([])
    else:
        pass
    plot_selected.set_ylabel(custom_ylabel)
    plot_selected.set_xlabel('Length [m]')
    """Print data"""
    # plot_to_compare.plot(length_steps_1, log_data_1[header], 'or', length_steps_2, log_data_2[header], 'ob', ms=0.2)
    plot_selected.plot(length_steps, log_data[headers])
    # plot_to_compare.plot(length_steps_1, log_data_1[header], 'or', length_steps_2, log_data_2[header], 'ob', ms=0.2)
    new_headers = list([])
    for word in headers:
        new_header = ''
        for char in word:
            if char == '>':
                char = ''
            new_header += char
        new_headers.append(new_header)
    plot_selected.legend(labels=new_headers, title='', loc='upper right',  borderaxespad=1, prop={'size': 6})
    plt.show()
    plt.show()


# РАЗМЕР ШАГА, ЭЛТЭЕ И КОЛ-ВО ШАГОВ + МАШИННОЕ ВРЕМЯ






#
# s1 = 'xxt_case1_v2_log'
# log_path_1 = os.path.join(os.getcwd(), 'log\{}.csv'.format(s1))
#
# s2 = 'xxt_case2_v1_log'
# log_path_2 = os.path.join(os.getcwd(), 'log\{}.csv'.format(s2))
#
# s3 = 'xxt_case3_v1_log'
# log_path_3 = os.path.join(os.getcwd(), 'log\{}.csv'.format(s3))
#
# s4 = 'xxt_case4_v1_log'
# log_path_4 = os.path.join(os.getcwd(), 'log\{}.csv'.format(s4))


# plot_compare(log_path_3, log_path_2, 'stepnum', ['1200 K', '1000K'])
# plot_compare(log_path_3, log_path_4, 'stepnum', ['1200 K', '1400K'], 'line')
# plot_compare(log_path_1, log_path_2, 'LTE', ['800 K', '1000K'])
# plot_compare(log_path_3, log_path_4, 'LTE', ['1200 K', '1400K'], 'points')

# plot_selected(log_path_3, ['step'])

# plot_profiles(log_path_3)





# plot_selected(log_path_3, ['C2H5* --> C2H4 + H*',
#                            'C2H6 + H* --> C2H5* + H2',
#                            'C2H4 + H* --> C2H5*',
#                            'C2H6 + CH3* --> C2H5* + CH4',
#                            'C2H6 --> CH3* + CH3*'])

# plot_selected(log_path_3, ['H*', 'CH3*', 'C2H5*'])
# plotlog('v22(big-fxd_1e-5)_log')
# plot_compare(log_path_1, log_path_2, ['C2H4'], ['Adaptive step', 'Fixed step (1e-5 m)'])
