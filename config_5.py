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
comp_x0 = dict({sp_db.comp0003.ID : 0.0176113723637944,
sp_db.comp0005.ID : 0.532084032272764,
sp_db.comp0007.ID : 0.00106463899490659,
sp_db.comp0006.ID : 6.01785422185644E-05,
sp_db.comp0011.ID : 0.00211998965131569,
sp_db.comp0009.ID : 0.00518430008977497,
sp_db.comp0013.ID : 0.003937357323084,
sp_db.comp0062.ID : 0.00101274699428162,
sp_db.comp0001.ID : 0.436925383767861
})


'''Setting up initial conditions'''
molflow = 489.5318442  # Feed stream molar flow [kgmol/hr]
P = 0.51 #* 0.101325  # Reaction Pressure [MPa]
T0 = 1050 + 273.15  # Initial Temperature [K] -previous 1050

'''Reactor rating'''
tube_L = 850  # Reaction tubes length [mm]
tube_ID = 24 # Reaction tubes Internal Diameter [mm]
tubes_No = 1  # Reaction tubes number

# pinter2004
