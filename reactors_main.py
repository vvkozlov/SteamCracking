'''
Header      : reactors_main.py
Created     : 22.06.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : 'Orchestrator program' for Plug-Flog Reactor simulation. Consolidate all scripts together.
'''


import os
import chemistry as rctr
import subprocess as sp
import reactors_config as cfg
import time

start_time = time.time()
print('Initializing calculations...')

'''Setting up reactions'''
rxnset = cfg.rxnset
rctr_compset = cfg.rctr_compset
comp_x0 = cfg.comp_x0

'''Setting up initial conditions'''
molflow = cfg.molflow  # Feed stream molar flow [kgmol/hr]
P = cfg.P  # Reaction Pressure [MPa]
T0 = cfg.T0  # Initial Temperature [K]

'''Reactor rating'''
tube_L = cfg.tube_L  # Reaction tubes length [mm]
tube_ID = cfg.tube_ID  # Reaction tubes Internal Diameter [mm]
tubes_No = cfg.tubes_No  # Reaction tubes number

'''Initializing feed Stream'''
inlet_stream = rctr.Stream(rctr_compset, comp_x0, molflow, P, T0)

'''Creating Reactor model'''
cstreactor = rctr.PFRreactor(tube_L / 1000, tube_ID / 1000, tubes_No, rxnset)
print('Starting calculations...')

'''Integrating through PFReactor model'''
outlet_stream, calc_hist = cstreactor.simulation(inlet_stream, 1e-2, True)
print('Calculations completed!')
print('runtime: {:.2f} ms'.format((time.time() - start_time) * 1000))
print('\nOutlet stream composition [mol. fract.]:\n\t', outlet_stream.COMPMOLFR)
print('Outlet stream temperature [K]:\n\t', outlet_stream.T)
print('Outlet stream act. vol. flow [m3/h]:\n\t', outlet_stream.FLVOLPR, '\n')

'''Saving results to .txt file'''
filename = 'log.txt'
filepath = os.path.join(os.getcwd(), filename)
print('Saving Results to {}...'.format(filepath))
log_file = open(filename, 'w')
log_file.write(calc_hist.to_string())
log_file.close()
sp.Popen(['notepad', os.path.join(os.getcwd(), filename)])
'''Saving results to .csv file'''
filename = 'log.csv'
calc_hist.to_csv(filename)

'''Plotting diagrams in matplotlib'''
print('Plotting graphs...')
rctr.plot_results(calc_hist)

'''Done)'''
print('done')