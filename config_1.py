'''
Header      : config_1.py
Created     : 08.08.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Contains all input data to simulate Plug-Flow Reactor. Added to get an opportunity to change
			  input quickly.
			  Input for base case with ethane and propane dehydrogenation reactions only.
'''


import databases.database_v1 as db


'''Setting up reagents reactions'''
PROPANE = db.PROPANE
PROPYLENE = db.PROPYLENE
ETHANE = db.ETHANE
ETHYLENE = db.ETHYLENE
H2 = db.H2

'''Setting up reactions themselves'''
rxn1 = db.rxn3
rxn2 = db.rxn4
rxnset = [rxn1, rxn2]

'''Compsets for individual rxns'''
rxn1_compset = rxn1.reagents
rxn2_compset = rxn2.reagents
'''Compset for all species in reactor'''
rctr_compset = list(set(rxn1_compset + rxn2_compset))

'''Feed stream mole composition [mol. frac.]'''
comp_x0 = dict({ETHANE.ID : 0.8, PROPANE.ID : 0.2, ETHYLENE.ID : 0, PROPYLENE.ID : 0, H2.ID: 0})

'''Setting up initial conditions'''
molflow = 1  # Feed stream molar flow [kgmol/hr]
P = 0.1  # Reaction Pressure [MPa]
T0 = 1000 + 273.15  # Initial Temperature [K]

'''Reactor rating'''
tube_L = 10000  # Reaction tubes length [mm]
tube_ID = 50  # Reaction tubes Internal Diameter [mm]
tubes_No = 1  # Reaction tubes number