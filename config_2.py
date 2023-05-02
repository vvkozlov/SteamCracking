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
comp_x0 = dict({sp_db.comp0005.ID : 1})

'''Setting up initial conditions'''
molflow = 9.034e-2  # Feed stream molar flow [kgmol/hr]
P = 3 * 0.101325  # Reaction Pressure [MPa]
T0 = 1000 + 273.15  # Initial Temperature [K]

'''Reactor rating'''
tube_L = 10000  # Reaction tubes length [mm]
tube_ID = 76  # Reaction tubes Internal Diameter [mm]
tubes_No = 1  # Reaction tubes number