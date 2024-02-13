'''
Header      : main.py
Created     : 22.06.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : 'Orchestrator program' for Plug-Flog Reactor simulation. Consolidate all scripts together.
'''
import os
import chemistry as rctr
from configs import config_1 as cfg
import plotter as pl
import time

'''Setting case  name'''
print('Save log and config files?: y/n')
save_option = str(input())
case_name = ''
if save_option == 'y':
      save_option = True
      print('Input case name:')
      case_name = str(input())
      '''Setting console clone file'''
      console_clone = pl.Logger(os.path.join(os.getcwd(), 'log\{}_status.dat'.format(case_name)))
else:
      save_option = False
      case_name = 'local case'

'''Logging execution details'''
start_time = time.time()
print('{}\tcase name: {}\t\tconfig file: {}'.format(time.strftime('%d.%m.%Y %H:%M', time.localtime(start_time + 14400)), case_name, cfg.__name__))  # + 14400 sec is conversion to GMT+7 time zone
print('-'*70, end= '\n\n')
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
print(f'Inlet stream temperature [K]:\t\t{inlet_stream.T : .3f}')
print(f'Inlet stream act. vol. flow [m3/hr]:\t{inlet_stream.FLVOL : .3f}')
print(f'Inlet stream molar flow [kgmol/hr]:\t{inlet_stream.FLMOL : .3f}')
print(f'Inlet stream mass flow [kg/hr]:\t\t{inlet_stream.FLMASS : .3f}\n')
inlet_mass = inlet_stream.FLMASS

'''Creating Reactor model'''
cstreactor = rctr.PFReactor(tube_L / 1000, tube_ID / 1000, tubes_No, rxnset)
print('Starting calculations...\n')

'''Setting denominator for number of outputed lines in logs'''
output_reduction = 1
'''Integrating through PFReactor model'''
outlet_stream, calc_hist = cstreactor.simulation(inlet_stream, 8e-5, True, output_reduction)
print('\nCalculations completed!')
'''Nice looking runtime report '''
runtime = time.time() - start_time
if runtime <= 60:
      print(f'runtime: {runtime % 60 : .3f} s\n')
elif runtime <=3600:
      print(f'runtime: {runtime // 60 : .0f} m {runtime % 60 : .1f} s\n')
elif runtime >3600:
      print(f'runtime: {runtime // 3600 : .0f} hr {(runtime - 3600 * runtime // 3600) // 60 : .0f} m {(runtime - 3600 * runtime // 3600) % 60 : .1f} s\n')
else:
      print(f'Waaaaaaay too long. Runtime: {runtime // 60 : .0f} m {runtime % 60 : .1f} s\n')

print(f'Outlet stream temperature [K]:\t\t{outlet_stream.T : .3f}')
print(f'Outlet stream act. vol. flow [m3/hr]:\t{outlet_stream.FLVOL : .3f}')
print(f'Outlet stream molar flow [m3/hr]:\t{outlet_stream.FLMOL : .3f}')
print(f'Outlet stream mass flow [kg/hr]:\t{outlet_stream.FLMASS : .3f}')
print(f'Stream mass flow divergence [kg/hr]:\t{(outlet_stream.FLMASS - inlet_mass) : .3f}\t'
      f'or [%] {(outlet_stream.FLMASS - inlet_mass) / inlet_mass * 100: .3f}\n')

if save_option:
      '''Saving results to .txt file'''
      log_filename = '{}_log'.format(case_name)
      filepath_txt = os.path.join(os.getcwd(), 'log\{}.txt'.format(log_filename))
      print('Saving Results to {}...'.format(filepath_txt))
      log_file = open(filepath_txt, 'w')
      log_file.write(calc_hist.to_string())
      log_file.close()
      # sp.Popen(['notepad', os.path.join(os.getcwd(), filename)])
      '''Saving results to .csv file'''
      filepath_csv = os.path.join(os.getcwd(), 'log\{}.csv'.format(log_filename))
      print('Saving Results to {}...'.format(filepath_csv))
      calc_hist.to_csv(filepath_csv)

      '''Plotting diagrams in matplotlib'''
      print('Plotting graphs...')
      # pl.plot_results(calc_hist)
      pl.plotlog(filepath_csv)
      console_clone.close()
else:
      pl.plothist(calc_hist)

'''Done)'''
print('done')
