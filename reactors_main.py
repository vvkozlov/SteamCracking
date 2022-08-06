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
rctr_compset = list(set(rxn1_compset + rxn2_compset))

'''Feed stream mole composition [mol. frac.]'''
comp_x0 = dict({ETHANE.name : 0.8, PROPANE.name : 0.2, ETHYLENE.name : 0, PROPYLENE.name : 0, H2.name: 0})

'''Setting up initial conditions'''
molflow = 1  # Feed stream molar flow [kgmol/hr]
P = 0.1  # Reaction Pressure [MPa]
T0 = 1000 + 273.15  # Initial Temperature [K]

'''Reactor rating'''
tube_L = 1000  # Reaction tubes length [mm]
tube_ID = 50  # Reaction tubes Internal Diameter [mm]
tubes_No = 1  # Reaction tubes number

'''Initializing feed Stream'''
inlet_stream = rctr.Stream(rctr_compset, comp_x0, molflow, P, T0)

'''Creating Reactor model'''
cstreactor = rctr.PFRreactor(tube_L / 1000, tube_ID / 1000, tubes_No, rxnset1)
print('Starting calculations...')

'''Integrating through PFReactor model'''
outlet_stream, calc_hist = cstreactor.simulation(inlet_stream, 1e-2, True)
print('Calculations completed!')

'''Saving results to .txt file'''
filename = 'log.txt'
filepath = os.path.join(os.getcwd(), filename)
print('Saving Results to {}...'.format(filepath))
log_file = open(filename, 'w')
log_file.write(calc_hist.to_string())
log_file.close()
sp.Popen(['notepad', os.path.join(os.getcwd(), filename)])

'''Plotting diagrams in matplotlib'''
print('Plotting graphs...')
rctr.plot_results(calc_hist)

'''Done)'''
print('done')
