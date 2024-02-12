"""
Header      : species_database_froment.py
Created     : 21.08.2023
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Storage of pure component's properties retrieved from one Soviet program
"""


from coreobjects import Species


"""
Example:
SPECIES = Species(ID, 'FORMULA', 'NAME', MW [kg/kgmol], [CPIG calc. option], [DCPIG eq. coeffs. x7], DHFORM [J/kgmol], PC [kPa], TC [degC], OMEGA)
"""

comp0002 = Species(2, 'H2', '', 2, 'Mayer-Kelly', [7.80000, -0.00249, 2.410E-06, -5.240E-10, 0, 0, 0, 0, 0, 0, 0], 0.000, 1357.755, -239.950, -0.220)
comp0003 = Species(3, 'CH4', '', 16, 'Mayer-Kelly', [0, 0.02520, -9.200E-06, 1.170E-09, 0, 0, 0, 0, 0, 0, 0], -74525022.200, 4640.685, -82.150, 0.008)
comp0004 = Species(4, 'C2H4', '', 28, 'Mayer-Kelly', [-2.55000, 0.04950, -3.500E-05, 1.040E-08, 0, 0, 0, 0, 0, 0, 0], 52334987.500, 5066.250, 8.850, 0.065)
comp0005 = Species(5, 'C2H6', '', 30, 'Mayer-Kelly', [-4.90000, 0.05980, -3.590E-05, 9.250E-09, 0, 0, 0, 0, 0, 0, 0], -84154659.900, 4883.865, 31.850, 0.098)
comp0006 = Species(6, 'C3H6', '', 42, 'Mayer-Kelly', [-4.35000, 0.07320, -4.800E-05, 1.320E-08, 0, 0, 0, 0, 0, 0, 0], 20305975.150, 4640.685, 90.850, 0.148)
comp0007 = Species(7, 'C3H8', '', 44, 'Mayer-Kelly', [-7.10000, 0.09420, -6.430E-05, 1.980E-08, 0, 0, 0, 0, 0, 0, 0], -102785915.450, 4255.650, 106.650, 0.152)
comp0008 = Species(8, 'C4H6', '', 54, 'Mayer-Kelly', [-1.38500, 0.08980, -7.090E-05, 2.300E-08, 0, 0, 0, 0, 0, 0, 0], 109694133.800, 4326.578, 151.850, 0.195)
comp0037 = Species(37, 'C2H2', '', 26, 'Mayer-Kelly', [9.15000, 0.01180, -7.200E-06, 2.390E-09, 0, 0, 0, 0, 0, 0, 0], 226087146.000, 6241.620, 35.850, 0.184)
comp0001 = Species(1, 'H2O', '', 18, 'Mayer-Kelly', [7.40000, 0.00160, 1.550E-06, -5.500E-10, 0, 0, 0, 0, 0, 0, 0], -241159622.400, 22088.850, 373.850, 0.344)


compset1 = [comp0001, comp0002,  comp0003,  comp0004,  comp0005,  comp0006,  comp0007,  comp0008,  comp0037]