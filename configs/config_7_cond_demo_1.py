'''
Header      :
Created     : 13.02.2024
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Kstovo 31-03-11 Regime  - "matches" with plant data @ 0.48 MW rctr heat duty

Reactor length set to 1 m to prepare set of demo cases (integration step vs. rctr conditions) for XXT
-> case 1 - 800 K
- case 2 - 1000 K
- case 3 - 1200 K
- case 4 - 1400 K
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
T0 = 800  # Initial Temperature [K]

'''Reactor rating'''
tube_L = 1000 # Reaction tubes length [mm]
tube_ID = 90 - 6.5 * 2 # Reaction tubes Internal Diameter [mm]
tubes_No = 1  # Reaction tubes number

'''Setting Furnace Heat Duty'''
furnace_duty = 0.48 * 10e3  # Furnace section heat duty [kJ/s]



