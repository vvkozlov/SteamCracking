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
import rctrs_engine_v3 as rctr
import database_v2 as db
import eos_engine_v4 as eos
import subprocess as sp
from eos_engine_v4 import flash_calc_PR_EOS_z_only as get_z_factor

'''Initializing Data'''
print('Initializing calculations...')
'''Setting up reagents for both reactions'''
PROPANE = db.PROPANE
PROPYLENE = db.PROPYLENE
ETHANE = db.ETHANE
ETHYLENE = db.ETHYLENE
H2 = db.H2

'''Setting up reactions themselves'''
rxn1 = db.rxn3
rxn2 = db.rxn4
rxnset1 = [rxn1, rxn2]

'''Compsets for individual rxns'''
rxn1_compset = rxn1.reagents
rxn2_compset = rxn2.reagents
'''Compset for all species in reactor'''
rctr_compset = set(rxn1_compset + rxn2_compset)



'''Feed stream mole composition [mol. frac.]'''
comp_x0 = dict({ETHANE.name : 0.8, PROPANE.name : 0.2, ETHYLENE.name : 0, PROPYLENE.name : 0, H2.name: 0})

'''Setting up initial conditions'''
molflow = 1  # Feed stream molar flow [kgmol/hr]
P = 0.1  # Reaction Pressure [MPa]
T0 = 1000 + 273.15  # Initial Temperature [K]

'''Resifence time for Plug-Flow Reactor application'''
tube_L = 1000  # Reaction tubes length [mm]
tube_ID = 50  # Reaction tubes Internal Diameter [mm]
tubes_No = 1  # Reaction tubes number

'''Residience time estimated with voluetric flowrate at inlet only'''
inlet_stream = rctr.Stream(rctr_compset, comp_x0, molflow, P, T0)
# print('Actual Volume Flow = {} [m3/h3]'.format(inlet_stream.FLVOLIG * get_z_factor(inlet_stream)))
print(inlet_stream.COMPMOLCIG)
print('Starting calculations...')
'''Creating Reactor model'''
cstreactor = rctr.PFRreactor(tube_L / 1000, tube_ID / 1000, tubes_No, rxnset1)

'''Integrating through PFReactor model'''
outlet_stream, calc_hist = cstreactor.Simulation(inlet_stream, 1e-6, True)
print('Calculations completed successfully!')

'''Saving results to .xlsx file'''
filepath = os.path.join(os.getcwd(), 'rctr_results.xlsx')
print('Saving Results to {}...'.format(filepath))
calc_hist.to_excel(filepath)
log = open('log.txt', 'w')
log.write(calc_hist.to_string())
log.close()
sp.Popen(['notepad', os.path.join(os.getcwd(), 'log.txt')])

'''Plotting diagrams in matplotlib'''
print('Plotting graphs...')
rctr.plot_results(calc_hist)

'''Done)'''
print('done')