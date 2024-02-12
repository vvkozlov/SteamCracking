"""
Header      : rxns_database_froment.py
Created     : 21.08.2023
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description :
"""




from databases.species_database_cowperthwaite import *
from chemistry import Reaction


"""
Example:
REACTION = Reaction('NAME', [Species xn], [stoic xn], [order [x], dH [kJ/mol], k0, E0 [kJ/mol], sequence)
If dH =0, rxn heat will be calculated from reagents enthalpy difference

sequence keys:
 1 - initiation rxn
 2 - propagation rxn
 3 - termination rxn
"""

rxn0001 = Reaction(1, 'C2H6 --> C2H4 + H2', [comp0005, comp0004, comp0002], [-1, 1, 1], [-1, 1, 1], 0, 4.6500E-10, 15.575,  1)
rxn0002 = Reaction(2, 'C2H4 + H2 --> C2H6', [comp0004, comp0002, comp0005], [-1, -1, 1], [-1, -1, 1], 0, 4.6500E-10, 15.575,  1)
rxn0003 = Reaction(3, 'C2H6 + C2H6 --> C3H8 + CH4', [comp0005, comp0005, comp0007, comp0003], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 3.8500E+08, 15.585,  1)
rxn0004 = Reaction(4, 'C3H6 --> C2H2 + CH4', [comp0006, comp0037, comp0003], [-1, 1, 1], [-1, 1, 1], 0, 9.8100E+05, 8.818,  1)
rxn0005 = Reaction(5, 'C2H2 + CH4 --> C3H6', [comp0037, comp0003, comp0006], [-1, -1, 1], [-1, -1, 1], 0, 9.8100E+05, 8.818,  2)
rxn0006 = Reaction(6, 'C2H2 + C2H4 --> C4H6', [comp0037, comp0004, comp0008], [-1, -1, 1], [-1, -1, 1], 0, 1.0300E+09, 9.855,  2)
rxn0007 = Reaction(7, 'C2H4 + C2H6 --> C3H6 + CH4', [comp0004, comp0005, comp0006, comp0003], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 7.0800E+10, 14.433,  2)


rxnset1 = [rxn0001, rxn0002, rxn0003, rxn0004, rxn0005, rxn0006, rxn0007]