"""
Header      : rxns_database_ussr.py
Created     : 26.04.2023
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

rxn0001 = Reaction(1, 'C2H6 --> CH3* + CH3*', [comp0005, comp0016, comp0016], [-1, 1, 1], [-1, 1, 1], 0, 4.0000E+16, 365.838,  1)
rxn0002 = Reaction(2, '1-C4H8 --> C3H5* + CH3*', [comp0009, comp0019, comp0016], [-1, 1, 1], [-1, 1, 1], 0, 8.0000E+16, 309.299,  1)
rxn0003 = Reaction(3, 'n-C4H10 --> C2H5* + C2H5*', [comp0011, comp0017, comp0017], [-1, 1, 1], [-1, 1, 1], 0, 1.5000E+16, 343.389,  1)
rxn0004 = Reaction(4, 'n-C4H10 --> n-C3H7* + CH3*', [comp0011, comp0020, comp0016], [-1, 1, 1], [-1, 1, 1], 0, 9.0000E+16, 357.524,  1)
rxn0005 = Reaction(5, 'C2H4 + H* --> C2H3* + H2', [comp0004, comp0015, comp0018, comp0002], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 8.0000E+11, 16.712,  2)
rxn0006 = Reaction(6, 'C2H6 + H* --> C2H5* + H2', [comp0005, comp0015, comp0017, comp0002], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 1.0000E+14, 40.575,  2)
rxn0007 = Reaction(7, 'C2H4 + CH3* --> C2H3* + CH4', [comp0004, comp0016, comp0018, comp0003], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 1.0000E+13, 54.377,  2)
rxn0008 = Reaction(8, 'C2H6 + CH3* --> C2H5* + CH4', [comp0005, comp0016, comp0017, comp0003], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 3.8000E+14, 69.010,  2)
rxn0009 = Reaction(9, 'C2H5* + H2 --> C2H6 + H*', [comp0017, comp0002, comp0005, comp0015], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 1.2000E+12, 58.617,  2)
rxn0010 = Reaction(10, 'C2H4 + C2H5* --> CH3* + C3H6', [comp0004, comp0017, comp0016, comp0006], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 3.0000E+12, 79.487,  2)
rxn0011 = Reaction(11, 'C2H2 + H* --> C2H3*', [comp0037, comp0015, comp0018], [-1, -1, 1], [-1, -1, 1], 0, 4.0000E+13, 5.438,  2)
rxn0012 = Reaction(12, 'C2H4 + H* --> C2H5*', [comp0004, comp0015, comp0017], [-1, -1, 1], [-1, -1, 1], 0, 1.0000E+13, 6.277,  2)
rxn0013 = Reaction(13, 'C3H6 + H* --> n-C3H7*', [comp0006, comp0015, comp0020], [-1, -1, 1], [-1, -1, 1], 0, 1.0000E+13, 12.139,  2)
rxn0014 = Reaction(14, 'C3H4 + H* --> C3H5*', [comp0036, comp0015, comp0019], [-1, -1, 1], [-1, -1, 1], 0, 3.5000E+13, 6.277,  2)
rxn0015 = Reaction(15, 'C4H6 + H* --> n-C4H7*', [comp0008, comp0015, comp0026], [-1, -1, 1], [-1, -1, 1], 0, 4.0000E+13, 5.438,  2)
rxn0016 = Reaction(16, 'C2H4 + CH3* --> n-C3H7*', [comp0004, comp0016, comp0020], [-1, -1, 1], [-1, -1, 1], 0, 2.0000E+11, 33.092,  2)
rxn0017 = Reaction(17, 'C2H4 + C2H3* --> n-C4H7*', [comp0004, comp0018, comp0026], [-1, -1, 1], [-1, -1, 1], 0, 5.0000E+10, 4.182,  2)
rxn0018 = Reaction(18, 'C4H6 + C2H3* --> CH3* + C5H6', [comp0008, comp0018, comp0016, comp0084], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 3.9800E+11, 20.869,  2)
rxn0019 = Reaction(19, 'C4H6 + C2H3* --> C6H6 + H2 + H*', [comp0008, comp0018, comp0050, comp0002, comp0015], [-1, -1, 1, 1, 1], [-1, -1, 1, 1, 1], 0, 1.2600E+11, 20.869,  2)
rxn0020 = Reaction(20, 'C2H4 + C2H5* --> n-C4H9*', [comp0004, comp0017, comp0022], [-1, -1, 1], [-1, -1, 1], 0, 1.5000E+10, 31.761,  2)
rxn0021 = Reaction(21, 'C3H6 + C2H5* --> n-C5H11*', [comp0006, comp0017, comp0041], [-1, -1, 1], [-1, -1, 1], 0, 1.3000E+10, 31.346,  2)
rxn0022 = Reaction(22, 'C2H4 + n-C3H7* --> n-C5H11*', [comp0004, comp0020, comp0041], [-1, -1, 1], [-1, -1, 1], 0, 2.0000E+10, 30.930,  2)
rxn0023 = Reaction(23, 'n-C4H9* --> sec-C4H9*', [comp0022, comp0023], [-1, 1], [-1, 1], 0, 3.1600E+12, 143.841,  2)
rxn0024 = Reaction(24, 'sec-C4H9* --> n-C4H9*', [comp0023, comp0022], [-1, 1], [-1, 1], 0, 3.1600E+12, 143.841,  2)
rxn0025 = Reaction(25, 'C2H3* --> C2H2 + H*', [comp0018, comp0037, comp0015], [-1, 1, 1], [-1, 1, 1], 0, 2.0000E+09, 132.201,  2)
rxn0026 = Reaction(26, 'C2H5* --> C2H4 + H*', [comp0017, comp0004, comp0015], [-1, 1, 1], [-1, 1, 1], 0, 3.2000E+13, 167.121,  2)
rxn0027 = Reaction(27, 'C3H5* --> C2H2 + CH3*', [comp0019, comp0037, comp0016], [-1, 1, 1], [-1, 1, 1], 0, 3.0000E+10, 151.324,  2)
rxn0028 = Reaction(28, 'C3H5* --> C3H4 + H*', [comp0019, comp0036, comp0015], [-1, 1, 1], [-1, 1, 1], 0, 8.9000E+12, 246.941,  2)
rxn0029 = Reaction(29, 'n-C3H7* --> C2H4 + CH3*', [comp0020, comp0004, comp0016], [-1, 1, 1], [-1, 1, 1], 0, 4.0000E+13, 136.358,  2)
rxn0030 = Reaction(30, 'n-C3H7* --> C3H6 + H*', [comp0020, comp0006, comp0015], [-1, 1, 1], [-1, 1, 1], 0, 2.0000E+13, 160.470,  2)
rxn0031 = Reaction(31, 'n-C4H7* --> C4H6 + H*', [comp0026, comp0008, comp0015], [-1, 1, 1], [-1, 1, 1], 0, 1.2000E+14, 232.806,  2)
rxn0032 = Reaction(32, 'n-C4H7* --> C2H4 + C2H3*', [comp0026, comp0004, comp0018], [-1, 1, 1], [-1, 1, 1], 0, 1.0000E+11, 133.032,  2)
rxn0033 = Reaction(33, 'n-C4H9* --> C2H4 + C2H5*', [comp0022, comp0004, comp0017], [-1, 1, 1], [-1, 1, 1], 0, 1.6000E+12, 116.403,  2)
rxn0034 = Reaction(34, 'n-C4H9* --> 1-C4H8 + H*', [comp0022, comp0009, comp0015], [-1, 1, 1], [-1, 1, 1], 0, 1.0000E+13, 149.661,  2)
rxn0035 = Reaction(35, 'sec-C4H9* --> C3H6 + CH3*', [comp0023, comp0006, comp0016], [-1, 1, 1], [-1, 1, 1], 0, 1.0000E+14, 141.347,  2)
rxn0036 = Reaction(36, 'n-C5H11* --> 1-C5H10 + H*', [comp0041, comp0053, comp0015], [-1, 1, 1], [-1, 1, 1], 0, 5.0000E+13, 149.661,  2)
rxn0037 = Reaction(37, 'n-C5H11* --> 1-C4H8 + CH3*', [comp0041, comp0009, comp0016], [-1, 1, 1], [-1, 1, 1], 0, 3.2000E+13, 124.718,  2)
rxn0038 = Reaction(38, 'n-C5H11* --> C2H4 + n-C3H7*', [comp0041, comp0004, comp0020], [-1, 1, 1], [-1, 1, 1], 0, 4.0000E+12, 116.403,  2)
rxn0039 = Reaction(39, 'C2H3* + H* --> C2H4', [comp0018, comp0015, comp0004], [-1, -1, 1], [-1, -1, 1], 0, 1.0000E+13, 0.000,  3)
rxn0040 = Reaction(40, 'C2H5* + H* --> C2H6', [comp0017, comp0015, comp0005], [-1, -1, 1], [-1, -1, 1], 0, 4.0000E+13, 0.000,  3)
rxn0041 = Reaction(41, 'C3H5* + H* --> C3H6', [comp0019, comp0015, comp0006], [-1, -1, 1], [-1, -1, 1], 0, 2.0000E+13, 0.000,  3)
rxn0042 = Reaction(42, 'n-C3H7* + H* --> C3H8', [comp0020, comp0015, comp0007], [-1, -1, 1], [-1, -1, 1], 0, 1.0000E+13, 0.000,  3)
rxn0043 = Reaction(43, 'n-C4H7* + H* --> 1-C4H8', [comp0026, comp0015, comp0009], [-1, -1, 1], [-1, -1, 1], 0, 2.0000E+13, 0.000,  3)
rxn0044 = Reaction(44, 'n-C4H9* + H* --> n-C4H10', [comp0022, comp0015, comp0011], [-1, -1, 1], [-1, -1, 1], 0, 1.0000E+13, 0.000,  3)
rxn0045 = Reaction(45, 'CH3* + CH3* --> C2H6', [comp0016, comp0016, comp0005], [-1, -1, 1], [-1, -1, 1], 0, 1.3000E+13, 0.000,  3)
rxn0046 = Reaction(46, 'C2H5* + CH3* --> C3H8', [comp0017, comp0016, comp0007], [-1, -1, 1], [-1, -1, 1], 0, 3.2000E+12, 0.000,  3)
rxn0047 = Reaction(47, 'C3H5* + CH3* --> 1-C4H8', [comp0019, comp0016, comp0009], [-1, -1, 1], [-1, -1, 1], 0, 3.2000E+12, 0.000,  3)
rxn0048 = Reaction(48, 'C2H3* + C2H3* --> C4H6', [comp0018, comp0018, comp0008], [-1, -1, 1], [-1, -1, 1], 0, 1.3000E+13, 0.000,  3)
rxn0049 = Reaction(49, 'C2H5* + C2H3* --> 1-C4H8', [comp0017, comp0018, comp0009], [-1, -1, 1], [-1, -1, 1], 0, 1.3000E+13, 0.000,  3)
rxn0050 = Reaction(50, 'C2H5* + C2H5* --> n-C4H10', [comp0017, comp0017, comp0011], [-1, -1, 1], [-1, -1, 1], 0, 4.0000E+11, 0.000,  3)
rxn0051 = Reaction(51, 'C2H5* + C2H5* --> C2H4 + C2H6', [comp0017, comp0017, comp0004, comp0005], [-1, -1, 1, 1], [-1, -1, 1, 1], 0, 5.0000E+10, 0.000,  3)
rxn0052 = Reaction(52, 'n-C4H7* + CH3* --> 1-C5H10', [comp0026, comp0016, comp0053], [-1, -1, 1], [-1, -1, 1], 0, 3.2000E+12, 0.000,  3)
rxn0053 = Reaction(53, 'n-C4H7* + C2H3* --> cyC6H10', [comp0026, comp0018, comp0078], [-1, -1, 1], [-1, -1, 1], 0, 1.3000E+13, 0.000,  3)
rxn0054 = Reaction(54, 'n-C4H7* + C2H5* --> 1-C6H12', [comp0026, comp0017, comp0062], [-1, -1, 1], [-1, -1, 1], 0, 3.2000E+12, 0.000,  3)
rxn0055 = Reaction(55, 'n-C4H7* + C3H5* --> n-C7H12', [comp0026, comp0019, comp0145], [-1, -1, 1], [-1, -1, 1], 0, 1.3000E+13, 0.000,  3)
rxn0056 = Reaction(56, 'n-C5H11* + H* --> n-C5H12', [comp0041, comp0015, comp0013], [-1, -1, 1], [-1, -1, 1], 0, 1.0000E+13, 0.000,  3)


rxnset1 = [rxn0001, rxn0002, rxn0003, rxn0004, rxn0005, rxn0006, rxn0007, rxn0008, rxn0009, rxn0010, rxn0011, rxn0012, rxn0013, rxn0014, rxn0015, rxn0016, rxn0017, rxn0018, rxn0019, rxn0020, rxn0021, rxn0022, rxn0023, rxn0024, rxn0025, rxn0026, rxn0027, rxn0028, rxn0029, rxn0030, rxn0031, rxn0032, rxn0033, rxn0034, rxn0035, rxn0036, rxn0037, rxn0038, rxn0039, rxn0040, rxn0041, rxn0042, rxn0043, rxn0044, rxn0045, rxn0046, rxn0047, rxn0048, rxn0049, rxn0050, rxn0051, rxn0052, rxn0053, rxn0054, rxn0055, rxn0056]