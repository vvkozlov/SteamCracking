'''
Header      : reactors_v5.py
Created     : 22.06.2022
Modified    : 24.07.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description :
Changes     :

'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import rctrs_engine as rctr
import database_v2 as db
import eos_engine as eos
import subprocess as sp
from eos_engine import flash_calc_PR_EOS_z_only as get_z_factor
import reactors_config as cfg

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
print('Outlet stream composition:\n', outlet_stream.COMPMOLFR)

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
