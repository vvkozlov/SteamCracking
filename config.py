'''
Header      : config.py
Created     : 08.08.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Contains all input data to simulate Plug-Flow Reactor. Added to get an opportunity to change
			  input quickly.
'''


import databases.species_database_ussr as sp_db
import databases.rxns_database_ussr as rxn_db


'''Compset for all species in reactor'''
rctr_compset = sp_db.compset1

'''Setting up reactions themselves'''
rxnset = rxn_db.rxnset1


'''Feed stream mole composition [mol. frac.]'''
comp_x0 = dict({sp_db.comp0005.name : 0.5, sp_db.comp0007.name : 0.5})

'''Setting up initial conditions'''
molflow = 1  # Feed stream molar flow [kgmol/hr]
P = 0.1  # Reaction Pressure [MPa]
T0 = 1000 + 273.15  # Initial Temperature [K]

'''Reactor rating'''
tube_L = 5000  # Reaction tubes length [mm]
tube_ID = 50  # Reaction tubes Internal Diameter [mm]
tubes_No = 1  # Reaction tubes number