'''
Header      : main.py
Created     : 22.06.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : 'Orchestrator program' for Plug-Flog Reactor simulation. Consolidate all scripts together.
'''
import math
import os
import chemistry as rctr
import subprocess as sp
import config_1 as cfg
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
print(f'Outlet stream temperature [K]:\t\t{inlet_stream.T : .3f}')
print(f'Inlet stream act. vol. flow [m3/hr]:\t{inlet_stream.FLVOL : .3f}')
print(f'Inlet stream mass flow [kg/hr]:\t\t{inlet_stream.FLMASS : .3f}\n')
inlet_mass = inlet_stream.FLMASS

'''Creating Reactor model'''
cstreactor = rctr.PFReactor(tube_L / 1000, tube_ID / 1000, tubes_No, rxnset)
# cstreactor.duty = 30
print('Starting calculations...\n')

'''Integrating through PFReactor model'''
outlet_stream, calc_hist = cstreactor.simulation(inlet_stream, 1e-3, True)
print('\nCalculations completed!')
'''Nice looking runtime report '''
runtime = time.time() - start_time
if runtime <= 60:
      print(f'runtime: {runtime % 60 : .3f} s\n')
elif runtime > 60:
      print(f'runtime: {runtime // 60 : .0f} m {runtime % 60 : .1f} s\n')
else:
      print(f'Waaaaaaay too long. Runtime: {runtime // 60 : .0f} m {runtime % 60 : .1f} s\n')

print(f'Outlet stream temperature [K]:\t\t{outlet_stream.T : .3f}')
print(f'Outlet stream act. vol. flow [m3/hr]:\t{outlet_stream.FLVOL : .3f}')
print(f'Outlet stream molar flow [m3/hr]:\t{outlet_stream.FLMOL : .3f}')
print(f'Outlet stream mass flow [kg/hr]:\t{outlet_stream.FLMASS : .3f}')
print(f'Stream mass divergence [kg/hr]:\t{(outlet_stream.FLMASS - inlet_mass) : .3f}\t'
      f'or [%] {(outlet_stream.FLMASS - inlet_mass) / inlet_mass * 100: .3f}\n')
# print('Outlet stream composition [mass. fract.]:\n\t', outlet_stream.COMPMASSFR)

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
# pl.plot_results(calc_hist)
pl.plotlog(filename)

'''Done)'''
print('done')
