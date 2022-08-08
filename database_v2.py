'''
Header      : database_v2.py
Created     : 02.08.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Storage of component's properties and reaction's parameters.
'''


from rctrs_engine import Species
from rctrs_engine import Reaction

NOCTANE = Species('n-C8H18', 114.23, [32.37317283, 105.8326168, 1635.6, 72.94353683, 746.4, 200, 1500], 0)
ISOOCTANE = Species('i-C8H18', 114.23, [0, 0, 0, 0, 0, 0, 0], 0)
NBUTANE = Species('C4H10', 58.12, [0, 0, 0, 0, 0, 0, 0], 0)
BUTENE1 = Species('C4H8', 56.11, [0, 0, 0, 0, 0, 0, 0], 0)
EBENZENE = Species('E-Benzene', 106.166000366211, [0, 0, 0, 0, 0, 0, 0], 0)
STYRENE = Species('Styrene', 104, [0, 0, 0, 0, 0, 0, 0], 0)
H2 = Species('Hydrogen', 2.01600003242493, [6.596207127, 2.283366772, 2466, 0.8980605713, 567.6, 250, 1500], 0)
BENZENE = Species('Benzene', 78.1100006103516, [0, 0, 0, 0, 0, 0, 0], 0)
ETHYLENE = Species('Ethylene', 28.0538005828857, [7.972676029, 22.64020254, 1596, 13.1604089, 740.8, 60, 1500], 0)
WATER = Species('Water', 18.0151004791260, [7.968615649, 6.398681571, 2610.5, 2.124773096, 1169, 100, 2273.15], 0)
METHANE = Species('Methane', 16.0429000854492, [0, 0, 0, 0, 0, 0, 0], 0)
CO2 = Species('CO2', 44.0097007751465, [0, 0, 0, 0, 0, 0, 0], 0)
CO = Species('CO', 28.0109004974365, [0, 0, 0, 0, 0, 0, 0], 0)
PROPYLENE = Species('Propene', 42.08064, [10.47387026, 35.97019203, 1398.8, 17.85468616, 616.46, 130, 1500], 0)
CL2 = Species('Cl2', 70.9054, [6.96044712, 2.191649947, 949, 2.395624343, 425, 50, 1500], 0)
C3H5CL = Species('AllylCl', 76.5254, [14.15615745, 34.3245438, 1677.7, 24.63217732, 769.64, 200, 1500], 0)
C3H6CL2 = Species('12-ClC3', 112.98604, [18.78714054, 41.62845132, 1715.7, 30.15907137, 765.1, 200, 1500], 0)
HCL = Species('HCl', 36.46064, [6.964029808, 2.161077673, 2093.8, -0.02555651094, 120, 50, 1500], 0)
ETHANE = Species('Ethane', 30.0699005126953, [10.570364, 20.23908474, 872.24, 16.03372504, 2430.4, 298.15, 1500], 0)
PROPANE = Species('Propane', 44.09652, [14.20512086, 30.24027897, 844.31, 20.58015668, 2482.7, 298.15, 1500], -104680000)
'''---------------------------------------------------------------------------'''
rxn1 = Reaction('Cl2 + Propene --> AllylCl + HCl', [CL2, PROPYLENE, C3H5CL, HCL], [-1, -1, 1, 1], -113.4919251, 1500000, 63.2672)
rxn2 = Reaction('Cl2 + Propene --> 12-ClC3', [CL2, PROPYLENE, C3H6CL2], [-1, -1, 1], -186.5439605, 90.46, 15.95636)
rxn3 = Reaction('C2H6 --> C2H4 + H2', [ETHANE, ETHYLENE, H2], [-1, 1, 1], 137.06780078125, 79432800000000, 297.2628)  # pre-exp is in [l/(mol*s)]
rxn4 = Reaction('C3H8 --> C2H6 + H2', [PROPANE, PROPYLENE, H2], [-1, 1, 1], 124.319, 5.01187E+13, 293.076)  # pre-exp is in [l/(mol*s)]