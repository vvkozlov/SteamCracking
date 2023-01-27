'''
Header      : main.py
Created     : 22.06.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : 'Orchestrator program' for Plug-Flog Reactor simulation. Consolidate all scripts together.
'''


import os
import chemistry as rctr
import subprocess as sp
import config as cfg
import plotter as pl
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
inlet_stream = rctr.Stream(rctr_compset, comp_x0, molflow, P, T0, 'IG')
print(f'Outlet stream mass flow [kg/h]:\t\t{inlet_stream.FLMASS : .3f}\n')

'''Creating Reactor model'''
cstreactor = rctr.PFReactor(tube_L / 1000, tube_ID / 1000, tubes_No, rxnset)
cstreactor.duty = 1
print('Starting calculations...\n')

'''Integrating through PFReactor model'''
outlet_stream, calc_hist = cstreactor.simulation(inlet_stream, 1e-2, True)
print('\nCalculations completed!')
print(f'runtime: {(time.time() - start_time) * 1000 : .2f} ms\n')
# print('Outlet stream composition [mol. fract.]:\n\t', outlet_stream.COMPMOLFR)
print(f'Outlet stream temperature [K]:\t\t{outlet_stream.T : .3f}')
print(f'Outlet stream act. vol. flow [m3/h]:\t{outlet_stream.FLVOL : .3f}')
print(f'Outlet stream mass flow [kg/h]:\t\t{outlet_stream.FLMASS : .3f}\n')

'''Saving results to .txt file'''
filename = 'log.txt'
filepath = os.path.join(os.getcwd(), filename)
print('Saving Results to {}...'.format(filepath))
log_file = open(filename, 'w')
log_file.write(calc_hist.to_string())
log_file.close()
# sp.Popen(['notepad', os.path.join(os.getcwd(), filename)])
'''Saving results to .csv file'''
filename = 'log.csv'
calc_hist.to_csv(filename)

'''Plotting diagrams in matplotlib'''
print('Plotting graphs...')
pl.plot_results(calc_hist)

'''Done)'''
print('done')