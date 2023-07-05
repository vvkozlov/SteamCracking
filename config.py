'''
Header      : config.py
Created     : 08.08.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Contains all input data to simulate Plug-Flow Reactor. Added to get an opportunity to change
			  input quickly.
'''


import databases.species_database_cowperthwaite as sp_db
import databases.rxns_database_cowperthwaite_v2 as rxn_db


'''Compset for all species in reactor'''
rctr_compset = sp_db.compset1

'''Setting up reactions themselves'''
rxnset = rxn_db.rxnset1


'''Feed stream mole composition [mol. frac.]'''
# comp_x0 = dict({sp_db.comp0003.ID : 0.4509, sp_db.comp0007.ID : 0.5491})
comp_x0 = dict({sp_db.comp0005.ID : 1})

'''Setting up initial conditions'''
# molflow = 9.034e-2  # Feed stream molar flow [kgmol/hr]
molflow = 100  # Feed stream molar flow [kgmol/hr]
# P = 7 * 0.101325  # Reaction Pressure [MPa]
P = 0.51 * 0.101325  # Reaction Pressure [MPa]
# T0 = 697 + 273.15  # Initial Temperature [K]
T0 = 1115 + 273.15  # Initial Temperature [K]

'''Reactor rating'''
# tube_L = 4000  # Reaction tubes length [mm]
tube_L = 4000  # Reaction tubes length [mm]
# tube_ID = 76  # Reaction tubes Internal Diameter [mm]
tube_ID = 76  # Reaction tubes Internal Diameter [mm]
# tubes_No = 1  # Reaction tubes number
tubes_No = 1  # Reaction tubes number

# Terrasug exmpl