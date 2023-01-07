'''
Header      : database_v2.py
Created     : 02.08.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Storage of component's properties and reaction's parameters.
'''
import pandas as pd

from rctrs_engine import Species
from rctrs_engine import Reaction

'''
Pure components properties database:
- DIPPR equation coefficients and DHFORM used from Aspen Plus V11
- PC, TC and OMEGA used from Aspen Hysys V10/V11 

Example:
SPECIES = Species(ID, 'NAME', 'FORMULA', MW [kg/kgmol], [CPIG calc. option], [DCPIG eq. coeffs. x7], DHFORM [J/kgmol], PC [kPa], TC [degC], OMEGA)
'''
NOCTANE = Species(-1, 'n-C8H18', '', 114.232002258301, 1, [32.37317283, 105.8326168, 1635.6, 72.94353683, 746.4, 200, 1500], 0, 2496.6201171875, 295.448022460938, 0.401800006628037)
ISOOCTANE = Species(-2, 'i-C8H18', '',  114.23, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
NBUTANE = Species(-3, 'C4H10', '',  58.1240005493164, 1, [0, 0, 0, 0, 0, 0, 0], 0, 3796.6201171875, 152.049005126953, 0.20100000500679)
BUTENE1 = Species(-4, 'C4H8', '',  56.11, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
EBENZENE = Species(-5, 'E-Benzene', '', 106.166000366211, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
STYRENE = Species(-6, 'Styrene', '', 104, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
H2 = Species(-7, 'Hydrogen', '', 2.01600003242493, 1, [6.596207127, 2.283366772, 2466, 0.8980605713, 567.6, 250, 1500], 0.0, 1315.5, -239.710001373291, -0.120090000331402)
BENZENE = Species(-8, 'Benzene', '', 78.1100006103516, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
ETHYLENE = Species(-9, 'Ethylene', '', 28.0538005828857, 1, [7.972676029, 22.64020254, 1596, 13.1604089, 740.8, 60, 1500], 52510000, 5031.79003906250, 9.20900878906252, 8.50000008940697e-002)
WATER = Species(-10, 'Water', '', 18.0151004791260, 1, [7.968615649, 6.398681571, 2610.5, 2.124773096, 1169, 100, 2273.15], 0, 0, 0, 0)
METHANE = Species(-11, 'Methane', '', 16.0429000854492, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
CO2 = Species(-12, 'CO2', '', 44.0097007751465, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
CO = Species(-13, 'CO', '', 28.0109004974365, 1, [0, 0, 0, 0, 0, 0, 0], 0, 0, 0, 0)
PROPYLENE = Species(-14, 'Propene', '', 42.08064, 1, [10.47387026, 35.97019203, 1398.8, 17.85468616, 616.46, 130, 1500], 20230000, 4620.41015625000, 91.85, 0.148000001907349)
CL2 = Species(-15, 'Cl2', '', 70.9054, 1, [6.96044712, 2.191649947, 949, 2.395624343, 425, 50, 1500], 0, 0, 0, 0)
C3H5CL = Species(-16, 'AllylCl', '', 76.5254, 1, [14.15615745, 34.3245438, 1677.7, 24.63217732, 769.64, 200, 1500], 0, 0, 0, 0)
C3H6CL2 = Species(-17, '12-ClC3', '', 112.98604, 1, [18.78714054, 41.62845132, 1715.7, 30.15907137, 765.1, 200, 1500], 0, 0, 0, 0)
HCL = Species(-18, 'HCl', '', 36.46064, 1, [6.964029808, 2.161077673, 2093.8, -0.02555651094, 120, 50, 1500], 0, 0, 0, 0)
ETHANE = Species(-19, 'Ethane', '', 30.0699005126953, 1, [10.570364, 20.23908474, 872.24, 16.03372504, 2430.4, 298.15, 1500], -83820000, 4883.85009765625, 32.2780090332031, 9.86000001430511e-002)
PROPANE = Species(-20, 'Propane', '', 44.09652, 1, [14.20512086, 30.24027897, 844.31, 20.58015668, 2482.7, 298.15, 1500], -104680000, 4256.66015625, 96.7480102539063, 0.152400001883507)

'''
Reactions data and properties database:

Example:
REACTION = Reaction('NAME', [Species xn], [stoic xn], dH [kJ/mol], k0, E0 [kJ/mol])
If dH =0, rxn heat will be calculated from reagents enthalpy difference
'''
rxn1 = Reaction('Cl2 + Propene --> AllylCl + HCl', [CL2, PROPYLENE, C3H5CL, HCL], [-1, -1, 1, 1], [-1, -1, 1, 1], -113.4919251, 1500000, 63.2672)
rxn2 = Reaction('Cl2 + Propene --> 12-ClC3', [CL2, PROPYLENE, C3H6CL2], [-1, -1, 1], [-1, -1, 1], -186.5439605, 90.46, 15.95636)
rxn3 = Reaction('C2H6 --> C2H4 + H2', [ETHANE, ETHYLENE, H2], [-1, 1, 1], [-1, 1, 1], 137.06780078125, 79432800000000, 297.2628)  # pre-exp is in [l/(mol*s)]
rxn4 = Reaction('C3H8 --> C2H6 + H2', [PROPANE, PROPYLENE, H2], [-1, 1, 1], [-1, 1, 1], 124.319, 5.01187E+13, 293.076)  # pre-exp is in [l/(mol*s)]
