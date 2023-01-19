"""
Header      : species_database_ussr.py
Created     : 08.01.2023
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Storage of pure component's properties retrieved from one Soviet program
"""


from coreobjects import Species


"""
Example:
SPECIES = Species(ID, 'NAME', 'FORMULA', MW [kg/kgmol], [CPIG calc. option], [DCPIG eq. coeffs. x7], DHFORM [J/kgmol], PC [kPa], TC [degC], OMEGA)
"""

comp0001 = Species(1, 'H2O', '', 18, 'Mayer-Kelly', [7.40000, 0.00160, 1.550E-06, -5.500E-10, 0, 0, 0, 0, 0, 0, 0], -241159622.400, 22088.850, 373.850, 0.344)
comp0002 = Species(2, 'H2', '', 2, 'Mayer-Kelly', [7.80000, -0.00249, 2.410E-06, -5.240E-10, 0, 0, 0, 0, 0, 0, 0], 0.000, 1357.755, -239.950, -0.220)
comp0003 = Species(3, 'CH4', '', 16, 'Mayer-Kelly', [0, 0.02520, -9.200E-06, 1.170E-09, 0, 0, 0, 0, 0, 0, 0], -74525022.200, 4640.685, -82.150, 0.008)
comp0004 = Species(4, 'C2H4', '', 28, 'Mayer-Kelly', [-2.55000, 0.04950, -3.500E-05, 1.040E-08, 0, 0, 0, 0, 0, 0, 0], 52334987.500, 5066.250, 8.850, 0.065)
comp0005 = Species(5, 'C2H6', '', 30, 'Mayer-Kelly', [-4.90000, 0.05980, -3.590E-05, 9.250E-09, 0, 0, 0, 0, 0, 0, 0], -84154659.900, 4883.865, 31.850, 0.098)
comp0006 = Species(6, 'C3H6', '', 42, 'Mayer-Kelly', [-4.35000, 0.07320, -4.800E-05, 1.320E-08, 0, 0, 0, 0, 0, 0, 0], 20305975.150, 4640.685, 90.850, 0.148)
comp0007 = Species(7, 'C3H8', '', 44, 'Mayer-Kelly', [-7.10000, 0.09420, -6.430E-05, 1.980E-08, 0, 0, 0, 0, 0, 0, 0], -102785915.450, 4255.650, 106.650, 0.152)
comp0008 = Species(8, 'C4H6', '', 54, 'Mayer-Kelly', [-1.38500, 0.08980, -7.090E-05, 2.300E-08, 0, 0, 0, 0, 0, 0, 0], 109694133.800, 4326.578, 151.850, 0.195)
comp0009 = Species(9, '1-C4H8', '', 56, 'Mayer-Kelly', [-5.18000, 0.10000, -6.730E-05, 1.910E-08, 0, 0, 0, 0, 0, 0, 0], -100231.968, 4022.603, 146.250, 0.187)
comp0010 = Species(10, 'i-C4H8', '', 56, 'Mayer-Kelly', [-1.96000, 0.09000, -5.600E-05, 1.490E-08, 0, 0, 0, 0, 0, 0, 0], -16914667.960, 4002.338, 148.550, 0.190)
comp0011 = Species(11, 'n-C4H10', '', 58, 'Mayer-Kelly', [-8.90000, 0.12200, -8.300E-05, 2.360E-08, 0, 0, 0, 0, 0, 0, 0], -126441329.800, 3799.688, 151.850, 0.193)
comp0012 = Species(12, 'i-C4H10', '', 58, 'Mayer-Kelly', [-8.90000, 0.12500, -8.750E-05, 2.560E-08, 0, 0, 0, 0, 0, 0, 0], -133977568.000, 3647.700, 138.750, 0.176)
comp0013 = Species(13, 'n-C5H12', '', 72, 'Mayer-Kelly', [-10.80000, 0.15100, -1.040E-04, 2.970E-08, 0, 0, 0, 0, 0, 0, 0], -146956644.900, 3374.123, 196.450, 0.251)
comp0014 = Species(14, 'i-C5H12', '', 72, 'Mayer-Kelly', [-11.85000, 0.15500, -1.090E-04, 3.180E-08, 0, 0, 0, 0, 0, 0, 0], -154492883.100, 3414.653, 187.650, 0.227)
comp0015 = Species(15, 'H*', '', 1, 'Mayer-Kelly', [7.80000, -0.00249, 2.410E-06, -5.240E-10, 0, 0, 0, 0, 0, 0, 0], 1534071615.000, 2968.820, 234.250, 0.296)
comp0016 = Species(16, 'CH3*', '', 15, 'Mayer-Kelly', [0, 0.02520, -9.200E-06, 1.170E-09, 0, 0, 0, 0, 0, 0, 0], -74525022.200, 4640.685, -82.150, 0.008)
comp0017 = Species(17, 'C2H5*', '', 29, 'Mayer-Kelly', [-4.90000, 0.05980, -3.590E-05, 9.250E-09, 0, 0, 0, 0, 0, 0, 0], -84154659.900, 4883.865, 31.850, 0.098)
comp0018 = Species(18, 'C2H3*', '', 27, 'Mayer-Kelly', [-2.55000, 0.04950, -3.500E-05, 1.040E-08, 0, 0, 0, 0, 0, 0, 0], 52334987.500, 5066.250, 8.850, 0.065)
comp0019 = Species(19, 'C3H5*', '', 41, 'Mayer-Kelly', [-4.35000, 0.07320, -4.800E-05, 1.320E-08, 0, 0, 0, 0, 0, 0, 0], 20305975.150, 4640.685, 90.850, 0.148)
comp0020 = Species(20, 'n-C3H7*', '', 43, 'Mayer-Kelly', [-7.10000, 0.09420, -6.430E-05, 1.980E-08, 0, 0, 0, 0, 0, 0, 0], -102785915.450, 4255.650, 106.650, 0.152)
comp0021 = Species(21, 'i-C3H7', '', 43, 'Mayer-Kelly', [-7.10000, 0.09420, -6.430E-05, 1.980E-08, 0, 0, 0, 0, 0, 0, 0], -102785915.450, 4255.650, 106.650, 0.152)
comp0022 = Species(22, 'n-C4H9*', '', 57, 'Mayer-Kelly', [-8.90000, 0.12200, -8.300E-05, 2.360E-08, 0, 0, 0, 0, 0, 0, 0], -126441329.800, 3799.688, 151.850, 0.193)
comp0023 = Species(23, 'sec-C4H9*', '', 57, 'Mayer-Kelly', [-8.90000, 0.12200, -8.300E-05, 2.360E-08, 0, 0, 0, 0, 0, 0, 0], -126441329.800, 3799.688, 151.850, 0.193)
comp0024 = Species(24, 'TP-C4H9*', '', 57, 'Mayer-Kelly', [-8.90000, 0.12200, -8.300E-05, 2.360E-08, 0, 0, 0, 0, 0, 0, 0], -126441329.800, 3799.688, 151.850, 0.193)
comp0025 = Species(25, 'i-C4H9*', '', 57, 'Mayer-Kelly', [-8.90000, 0.12500, -8.750E-05, 2.560E-08, 0, 0, 0, 0, 0, 0, 0], -133977568.000, 3647.700, 138.750, 0.176)
comp0026 = Species(26, 'n-C4H7*', '', 55, 'Mayer-Kelly', [-5.18000, 0.10000, -6.730E-05, 1.910E-08, 0, 0, 0, 0, 0, 0, 0], -100231.968, 4022.603, 146.250, 0.187)
comp0027 = Species(27, 'OH*', '', 17, 'Mayer-Kelly', [7.40000, 0.00160, 1.550E-06, -5.500E-10, 0, 0, 0, 0, 0, 0, 0], -143311300.800, 2968.820, 234.250, 0.296)
comp0028 = Species(28, 'O**', '', 16, 'Mayer-Kelly', [6.71300, 0.00000, 4.170E-06, -2.540E-09, 0, 0, 0, 0, 0, 0, 0], 0.000, 5045.985, -118.550, 0.021)
comp0029 = Species(29, 'CH3O*', '', 31, 'Mayer-Kelly', [1.84340, 0.00353, -2.400E-05, 5.685E-09, 0, 0, 0, 0, 0, 0, 0], -166467128.240, 5572.875, 187.850, 0.303)
comp0030 = Species(30, 'HCO*', '', 29, 'Mayer-Kelly', [1.39000, 0.02140, -9.600E-06, 1.360E-09, 0, 0, 0, 0, 0, 0, 0], -108438094.100, 4053.000, 226.850, 0.253)
comp0031 = Species(31, 'OCOH*', '', 45, 'Mayer-Kelly', [1.52530, 0.00913, -2.590E-05, 6.445E-09, 0, 0, 0, 0, 0, 0, 0], -61127265.400, 6484.800, 106.850, 0.207)
comp0032 = Species(32, 'H2O2', '', 34, 'Mayer-Kelly', [4.25000, 0.02690, -2.500E-05, 8.850E-09, 0, 0, 0, 0, 0, 0, 0], -136070967.500, 15198.750, 326.850, 0.090)
comp0033 = Species(33, 'H-CH=O', '', 30, 'Mayer-Kelly', [1.39000, 0.02140, -9.600E-06, 1.360E-09, 0, 0, 0, 0, 0, 0, 0], -108438094.100, 4053.000, 226.850, 0.253)
comp0034 = Species(34, 'CO', '', 28, 'Mayer-Kelly', [6.70000, -0.00079, -3.850E-06, 1.860E-09, 0, 0, 0, 0, 0, 0, 0], -110531493.600, 3495.713, -140.450, 0.049)
comp0035 = Species(35, 'CO2', '', 44, 'Mayer-Kelly', [5.62000, 0.01430, -9.400E-06, 2.440E-09, 0, 0, 0, 0, 0, 0, 0], -393559106.000, 7386.593, 27.850, 0.225)
comp0036 = Species(36, 'C3H4', 'Methylacetylene', 40, 'Mayer-Kelly', [3.96000, 0.04280, -2.600E-05, 6.880E-09, 0, 0, 0, 0, 0, 0, 0], 184219156.000, 5349.960, 127.850, 0.218)
comp0037 = Species(37, 'C2H2', '', 26, 'Mayer-Kelly', [9.15000, 0.01180, -7.200E-06, 2.390E-09, 0, 0, 0, 0, 0, 0, 0], 226087146.000, 6241.620, 35.850, 0.184)
comp0038 = Species(38, '1,7-C8H14', '', 110, 'Mayer-Kelly', [9.00000, 0.12600, -4.560E-05, 0, 0, 0, 0, 0, 0, 0, 0], 39146570.650, 2674.980, 298.850, 0.400)
comp0039 = Species(39, '2-C4H8', '', 56, 'Mayer-Kelly', [-5.20000, 0.09600, -6.120E-05, 1.660E-08, 0, 0, 0, 0, 0, 0, 0], -11178753.330, 4154.325, 156.850, 0.210)
comp0040 = Species(40, 'i-C4H7*', '', 55, 'Mayer-Kelly', [-1.96000, 0.09000, -5.600E-05, 1.490E-08, 0, 0, 0, 0, 0, 0, 0], -16914667.960, 4002.338, 148.550, 0.190)
comp0041 = Species(41, 'n-C5H11*', '', 71, 'Mayer-Kelly', [-10.80000, 0.15100, -1.040E-04, 2.970E-08, 0, 0, 0, 0, 0, 0, 0], -146956644.900, 3374.123, 196.450, 0.251)
comp0042 = Species(42, 'i-C5H11*', '', 71, 'Mayer-Kelly', [-11.85000, 0.15500, -1.090E-04, 3.180E-08, 0, 0, 0, 0, 0, 0, 0], -154492883.100, 3414.653, 187.650, 0.227)
comp0043 = Species(43, 'C6H12**', '', 84, 'Mayer-Kelly', [-8.75000, 0.16600, -1.060E-04, 2.960E-08, 0, 0, 0, 0, 0, 0, 0], -41867990.000, 3242.400, 230.850, 0.285)
comp0044 = Species(44, '2-C6H11*', '', 83, 'Mayer-Kelly', [-16.30000, 0.17800, -1.290E-04, 3.900E-08, 0, 0, 0, 0, 0, 0, 0], -46473468.900, 3394.388, 243.850, 0.250)
comp0045 = Species(45, 'cyC6H11*', '', 83, 'Mayer-Kelly', [-35.00000, 0.21500, -1.570E-04, 4.300E-08, 0, 0, 0, 0, 0, 0, 0], -123510570.500, 4053.000, 279.850, 0.213)
comp0046 = Species(46, 'cyC6H9*', '', 81, 'Mayer-Kelly', [-21.30000, 0.18700, -1.460E-04, 4.400E-08, 0, 0, 0, 0, 0, 0, 0], -5400970.710, 4235.385, 285.850, 0.210)
comp0047 = Species(47, 'cyC6H7*', '', 79, 'Mayer-Kelly', [3.54000, 0.00441, 4.510E-06, 0, 0, 0, 0, 0, 0, 0, 0], 112876101.040, 3546.375, 259.850, 0.317)
comp0048 = Species(48, 'C4H5*', '', 53, 'Mayer-Kelly', [-1.38500, 0.08980, -7.090E-05, 2.300E-08, 0, 0, 0, 0, 0, 0, 0], 109694133.800, 4326.578, 151.850, 0.195)
comp0049 = Species(49, '1,3-C5H8', '', 68, 'Mayer-Kelly', [-7.88000, 0.12800, -9.820E-05, 3.100E-08, 0, 0, 0, 0, 0, 0, 0], 78293141.300, 3850.350, 214.850, 0.164)
comp0050 = Species(50, 'C6H6', '', 78, 'Mayer-Kelly', [-14.60000, 0.14000, -1.080E-04, 3.300E-08, 0, 0, 0, 0, 0, 0, 0], 82479940.300, 4924.395, 288.850, 0.212)
comp0051 = Species(51, 'cyC6H12', '', 84, 'Mayer-Kelly', [-35.00000, 0.21500, -1.570E-04, 4.300E-08, 0, 0, 0, 0, 0, 0, 0], -123510570.500, 4053.000, 279.850, 0.213)
comp0052 = Species(52, 'cyC6H8', '1,3-Cyclohexadiene', 80, 'Mayer-Kelly', [3.54000, 0.00441, 4.510E-06, 0, 0, 0, 0, 0, 0, 0, 0], 112876101.040, 3546.375, 259.850, 0.317)
comp0053 = Species(53, '1-C5H10', '', 70, 'Mayer-Kelly', [-6.82000, 0.13000, -8.500E-05, 2.360E-08, 0, 0, 0, 0, 0, 0, 0], -20933995.000, 4053.000, 191.550, 0.245)
comp0054 = Species(54, '2-C5H10', '', 70, 'Mayer-Kelly', [-13.40000, 0.14600, -1.050E-04, 3.150E-08, 0, 0, 0, 0, 0, 0, 0], -30144952.800, 4093.530, 202.350, 0.240)
comp0055 = Species(55, '1-C5H9*', '', 69, 'Mayer-Kelly', [-6.82000, 0.13000, -8.500E-05, 2.360E-08, 0, 0, 0, 0, 0, 0, 0], -20933995.000, 4053.000, 191.550, 0.245)
comp0056 = Species(56, '2-C5H9*', '', 69, 'Mayer-Kelly', [-13.40000, 0.14600, -1.050E-04, 3.150E-08, 0, 0, 0, 0, 0, 0, 0], -30144952.800, 4093.530, 202.350, 0.240)
comp0057 = Species(57, '3,3-C7H15*', '', 99, 'Mayer-Kelly', [-13.75000, 0.20500, -1.410E-04, 4.030E-08, 0, 0, 0, 0, 0, 0, 0], -201803711.800, 3647.700, 262.850, 0.270)
comp0058 = Species(58, 'n-C6H13', '', 85, 'Mayer-Kelly', [-13.05000, 0.18200, -1.260E-04, 3.650E-08, 0, 0, 0, 0, 0, 0, 0], -167471960.000, 3029.618, 234.650, 0.296)
comp0059 = Species(59, '2-C6H13*', '', 85, 'Mayer-Kelly', [-14.80000, 0.19300, -1.420E-04, 4.320E-08, 0, 0, 0, 0, 0, 0, 0], -175845558.000, 3029.618, 224.650, 0.279)
comp0060 = Species(60, '3-C6H13*', '', 85, 'Mayer-Kelly', [-13.70000, 0.18400, -1.290E-04, 3.770E-08, 0, 0, 0, 0, 0, 0, 0], -173752158.500, 3120.810, 231.450, 0.275)
comp0061 = Species(61, 'n-C7H16', '', 100, 'Mayer-Kelly', [-13.90000, 0.20600, -1.420E-04, 4.050E-08, 0, 0, 0, 0, 0, 0, 0], -188405955.000, 2735.775, 266.950, 0.351)
comp0062 = Species(62, '1-C6H12', '', 84, 'Mayer-Kelly', [-8.75000, 0.16600, -1.060E-04, 2.960E-08, 0, 0, 0, 0, 0, 0, 0], -41867990.000, 3242.400, 230.850, 0.285)
comp0063 = Species(63, 'n-C7H15*', '', 99, 'Mayer-Kelly', [-13.90000, 0.20600, -1.420E-04, 4.050E-08, 0, 0, 0, 0, 0, 0, 0], -188405955.000, 2735.775, 266.950, 0.351)
comp0064 = Species(64, '1-C6H11*', '', 83, 'Mayer-Kelly', [-8.75000, 0.16600, -1.060E-04, 2.960E-08, 0, 0, 0, 0, 0, 0, 0], -41867990.000, 3242.400, 230.850, 0.285)
comp0065 = Species(65, '2-C7H15*', '', 99, 'Mayer-Kelly', [-13.90000, 0.20600, -1.420E-04, 4.050E-08, 0, 0, 0, 0, 0, 0, 0], -195104833.400, 2756.040, 257.850, 0.330)
comp0066 = Species(66, '2-C6H12', '', 84, 'Mayer-Kelly', [-16.30000, 0.17800, -1.290E-04, 3.900E-08, 0, 0, 0, 0, 0, 0, 0], -46473468.900, 3394.388, 243.850, 0.250)
comp0067 = Species(67, 'i-C5H10', '', 70, 'Mayer-Kelly', [6.10000, 0.09400, -5.120E-05, 1.150E-08, 0, 0, 0, 0, 0, 0, 0], -29181989.030, 3434.918, 190.850, 0.260)
comp0068 = Species(68, '3-C7H16', '', 100, 'Mayer-Kelly', [-13.90000, 0.20600, -1.420E-04, 4.050E-08, 0, 0, 0, 0, 0, 0, 0], -193011433.900, 2756.040, 257.850, 0.324)
comp0069 = Species(69, '3-C7H15*', '', 99, 'Mayer-Kelly', [-13.90000, 0.20600, -1.420E-04, 4.050E-08, 0, 0, 0, 0, 0, 0, 0], -193011433.900, 2756.040, 257.850, 0.324)
comp0070 = Species(70, 'i-C6H12', '', 84, 'Mayer-Kelly', [-21.60000, 0.20400, -1.660E-04, 5.480E-08, 0, 0, 0, 0, 0, 0, 0], -50660267.900, 3394.388, 243.850, 0.260)
comp0071 = Species(71, 'n-C6H14', '', 86, 'Mayer-Kelly', [-13.05000, 0.18200, -1.260E-04, 3.650E-08, 0, 0, 0, 0, 0, 0, 0], -167471960.000, 3029.618, 234.650, 0.296)
comp0072 = Species(72, '2-C6H14', '', 86, 'Mayer-Kelly', [-14.80000, 0.19300, -1.420E-04, 4.320E-08, 0, 0, 0, 0, 0, 0, 0], -175845558.000, 3029.618, 224.650, 0.279)
comp0073 = Species(73, '3-C6H14', '', 86, 'Mayer-Kelly', [-13.70000, 0.18400, -1.290E-04, 3.770E-08, 0, 0, 0, 0, 0, 0, 0], -173752158.500, 3120.810, 231.450, 0.275)
comp0074 = Species(74, '2-C7H16', '', 100, 'Mayer-Kelly', [-13.90000, 0.20600, -1.420E-04, 4.050E-08, 0, 0, 0, 0, 0, 0, 0], -195104833.400, 2756.040, 257.850, 0.330)
comp0075 = Species(75, '1,3-C6H10', '', 82, 'Mayer-Kelly', [5.70000, 0.08640, -2.760E-05, 0, 0, 0, 0, 0, 0, 0, 0], 74525022.200, 3293.063, 229.850, 0.160)
comp0076 = Species(76, '3,3-C7H16', '', 100, 'Mayer-Kelly', [-13.75000, 0.20500, -1.410E-04, 4.030E-08, 0, 0, 0, 0, 0, 0, 0], -201803711.800, 3647.700, 262.850, 0.270)
comp0077 = Species(77, 'i-C5H9*', '', 69, 'Mayer-Kelly', [6.10000, 0.09400, -5.120E-05, 1.150E-08, 0, 0, 0, 0, 0, 0, 0], -29181989.030, 3434.918, 190.850, 0.260)
comp0078 = Species(78, 'cyC6H10', '', 82, 'Mayer-Kelly', [-21.30000, 0.18700, -1.460E-04, 4.400E-08, 0, 0, 0, 0, 0, 0, 0], -5400970.710, 4235.385, 285.850, 0.210)
comp0079 = Species(79, 'C7H14', 'Methylcyclohexane', 98, 'Mayer-Kelly', [-31.80000, 0.24500, -1.730E-04, 4.820E-08, 0, 0, 0, 0, 0, 0, 0], -154911563.000, 3475.448, 298.850, 0.358)
comp0080 = Species(80, 'C7H8', '', 92, 'Mayer-Kelly', [-18.45000, 0.17300, -1.310E-04, 3.950E-08, 0, 0, 0, 0, 0, 0, 0], 50032248.050, 4215.120, 318.550, 0.257)
comp0081 = Species(81, 'C8H10', '', 106, 'Mayer-Kelly', [-17.50000, 0.19400, -1.460E-04, 4.400E-08, 0, 0, 0, 0, 0, 0, 0], 29726272.900, 3850.350, 346.350, 0.301)
comp0082 = Species(82, 'C8H8', 'Styrene', 104, 'Mayer-Kelly', [-18.10000, 0.25300, -1.930E-04, 3.990E-08, 0, 0, 0, 0, 0, 0, 0], 147375324.800, 4053.000, 272.850, 0.257)
comp0083 = Species(83, 'biC14H14', 'Bybenzyl', 182, 'Mayer-Kelly', [4.82000, 0.17800, -6.020E-05, 0, 0, 0, 0, 0, 0, 0, 0], 51497627.700, 4154.325, 541.850, 0.375)
comp0084 = Species(84, 'C5H6', 'Cyclopentadiene', 66, 'Mayer-Kelly', [3.42000, 0.04890, -4.960E-06, 0, 0, 0, 0, 0, 0, 0, 0], 133140208.200, 4559.625, 236.850, 0.132)
comp0085 = Species(85, 'cyC6H8', 'Methylcyclopentadiene', 80, 'Mayer-Kelly', [3.58000, 0.06470, -2.320E-05, 0, 0, 0, 0, 0, 0, 0, 0], 113880932.800, 4255.650, 266.850, 0.223)
comp0086 = Species(86, 'cyC7H10', '1-Methylcyclohexadiene', 94, 'Mayer-Kelly', [2.31000, 0.08690, -1.650E-05, 0, 0, 0, 0, 0, 0, 0, 0], 80386540.800, 4053.000, 286.850, 0.492)
comp0087 = Species(87, 'dicyC10H10', '2,3-Dihydronaphthalene', 130, 'Mayer-Kelly', [3.10000, 0.10900, -3.420E-05, 0, 0, 0, 0, 0, 0, 0, 0], 123929250.400, 3212.003, 451.850, 0.292)
comp0088 = Species(88, 'i-C6H11*', '', 83, 'Mayer-Kelly', [-21.60000, 0.20400, -1.660E-04, 5.480E-08, 0, 0, 0, 0, 0, 0, 0], -50660267.900, 3394.388, 243.850, 0.260)
comp0089 = Species(89, 'C5H10**', '', 70, 'Mayer-Kelly', [-6.82000, 0.13000, -8.500E-05, 2.360E-08, 0, 0, 0, 0, 0, 0, 0], -20933995.000, 4053.000, 191.550, 0.245)
comp0090 = Species(90, 'C7H14**', '', 98, 'Mayer-Kelly', [-31.80000, 0.24500, -1.730E-04, 4.820E-08, 0, 0, 0, 0, 0, 0, 0], -154911563.000, 3475.448, 298.850, 0.358)
comp0091 = Species(91, 'C7H13*', '', 97, 'Mayer-Kelly', [-31.80000, 0.24500, -1.730E-04, 4.820E-08, 0, 0, 0, 0, 0, 0, 0], -154911563.000, 3475.448, 298.850, 0.358)
comp0092 = Species(92, 'C7H7*', '', 91, 'Mayer-Kelly', [-18.45000, 0.17300, -1.310E-04, 3.950E-08, 0, 0, 0, 0, 0, 0, 0], 50032248.050, 4215.120, 318.550, 0.257)
comp0093 = Species(93, 'C8H9*', '', 105, 'Mayer-Kelly', [-17.50000, 0.19400, -1.460E-04, 4.400E-08, 0, 0, 0, 0, 0, 0, 0], 29726272.900, 3850.350, 346.350, 0.301)
comp0094 = Species(94, 'C7H9*', '', 93, 'Mayer-Kelly', [2.31000, 0.08690, -1.650E-05, 0, 0, 0, 0, 0, 0, 0, 0], 80386540.800, 4053.000, 286.850, 0.492)
comp0095 = Species(95, 'C5H8', '', 68, 'Mayer-Kelly', [3.02000, 0.06780, -1.480E-05, 0, 0, 0, 0, 0, 0, 0, 0], 39774590.500, 4053.000, 226.850, 0.193)
comp0096 = Species(96, 'cyC7H12', 'Methylcyclohexene', 96, 'Mayer-Kelly', [-14.88000, 0.19300, -1.380E-04, 0, 0, 0, 0, 0, 0, 0, 0], -35964603.410, 4053.000, 286.850, 0.490)
comp0097 = Species(97, 'C6H9*', '', 81, 'Mayer-Kelly', [5.70000, 0.08640, -2.760E-05, 0, 0, 0, 0, 0, 0, 0, 0], 74525022.200, 3293.063, 229.850, 0.160)
comp0098 = Species(98, 'C5H7*', '', 67, 'Mayer-Kelly', [3.02000, 0.06780, -1.480E-05, 0, 0, 0, 0, 0, 0, 0, 0], 39774590.500, 4053.000, 226.850, 0.193)
comp0099 = Species(99, 'C7H11*', '', 95, 'Mayer-Kelly', [-14.88000, 0.19300, -1.380E-04, 0, 0, 0, 0, 0, 0, 0, 0], -35964603.410, 4053.000, 286.850, 0.490)
comp0100 = Species(100, '3-C7H14', '', 98, 'Mayer-Kelly', [6.82000, 0.12500, -4.610E-05, 0, 0, 0, 0, 0, 0, 0, 0], -62341437.110, 3475.448, 298.950, 0.182)
comp0101 = Species(101, '3-C7H12', '3-Methylhexadiene', 96, 'Mayer-Kelly', [7.45000, 0.10200, -3.300E-05, 0, 0, 0, 0, 0, 0, 0, 0], 24702114.100, 3445.050, 304.850, 0.210)
comp0102 = Species(102, 'n-C10H22', '', 142, 'Mayer-Kelly', [8.18000, 0.19400, -7.050E-05, 0, 0, 0, 0, 0, 0, 0, 0], -249826296.330, 2107.560, 344.450, 0.490)
comp0103 = Species(103, 'n-C10H20', '', 140, 'Mayer-Kelly', [9.74000, 0.17900, -6.500E-05, 0, 0, 0, 0, 0, 0, 0, 0], -124222326.330, 2208.885, 341.850, 0.491)
comp0104 = Species(104, 'i-C10H20', '', 140, 'Mayer-Kelly', [10.15000, 0.18100, -6.670E-05, 0, 0, 0, 0, 0, 0, 0, 0], -136489647.400, 2127.825, 339.850, 0.533)
comp0105 = Species(105, 'n-C13H28', '', 184, 'Mayer-Kelly', [11.51000, 0.24900, -9.110E-05, 0, 0, 0, 0, 0, 0, 0, 0], -320708803.400, 1722.525, 384.650, 0.623)
comp0106 = Species(106, '5-C14H30', '', 198, 'Mayer-Kelly', [13.03000, 0.26900, -9.970E-05, 0, 0, 0, 0, 0, 0, 0, 0], -349179036.600, 1621.200, 420.850, 0.630)
comp0107 = Species(107, '1-C14H28', '', 196, 'Mayer-Kelly', [14.18000, 0.25300, -9.240E-05, 0, 0, 0, 0, 0, 0, 0, 0], -206660398.640, 1560.405, 415.850, 0.644)
comp0108 = Species(108, 'n-C17H36', '', 240, 'Mayer-Kelly', [15.95000, 0.32300, -1.190E-04, 0, 0, 0, 0, 0, 0, 0, 0], -407794222.600, 1317.225, 459.850, 0.770)
comp0109 = Species(109, 'i-C17H36', '', 240, 'Mayer-Kelly', [16.36000, 0.32500, -1.200E-04, 0, 0, 0, 0, 0, 0, 0, 0], -414493101.000, 1317.225, 459.850, 0.760)
comp0110 = Species(110, 'C15H30', 'Cyclopentadecane', 210, 'Mayer-Kelly', [12.58000, 0.27500, -9.970E-05, 0, 0, 0, 0, 0, 0, 0, 0], -333478540.350, 1459.080, 426.850, 0.660)
comp0111 = Species(111, 'cyC8H16', 'Dimethylcyclohexane', 112, 'Mayer-Kelly', [4.62000, 0.14200, -4.770E-05, 0, 0, 0, 0, 0, 0, 0, 0], -173333478.600, 3039.750, 326.850, 0.235)
comp0112 = Species(112, 'C10H18', 'Cyclodecane', 138, 'Mayer-Kelly', [6.73000, 0.18000, -6.250E-05, 0, 0, 0, 0, 0, 0, 0, 0], -218132227.900, 3039.750, 396.850, 0.295)
comp0113 = Species(113, 'C15H22', 'TETPЛ (?)', 202, 'Mayer-Kelly', [-10.34000, 0.29500, -1.170E-04, 0, 0, 0, 0, 0, 0, 0, 0], -119868055.370, 3039.750, 566.850, 0.380)
comp0114 = Species(114, 'C14H18', 'TETPAЛ (?)', 186, 'Mayer-Kelly', [-7.85000, 0.26200, -1.100E-04, 0, 0, 0, 0, 0, 0, 0, 0], 37388115.070, 3546.375, 546.850, 0.440)
comp0115 = Species(115, 'C12H14', '1,4-Dimethylnaphthalene', 158, 'Mayer-Kelly', [-10.34000, 0.22300, -8.890E-05, 0, 0, 0, 0, 0, 0, 0, 0], 95877697.100, 3039.750, 516.850, 0.590)
comp0116 = Species(116, 'C9H10', 'alpha-Methylstyrene', 118, 'Mayer-Kelly', [-7.45000, 0.10400, -3.150E-05, 0, 0, 0, 0, 0, 0, 0, 0], 116393012.200, 3404.520, 390.850, 0.244)
comp0117 = Species(117, 'C10H8', 'Naphthalene', 128, 'Mayer-Kelly', [3.15000, 0.10900, -3.480E-05, 0, 0, 0, 0, 0, 0, 0, 0], 151980803.700, 3982.073, 475.250, 0.302)
comp0118 = Species(118, 'C12H12', '1,5-Dimethylnaphthalene', 156, 'Mayer-Kelly', [4.79000, 0.14200, -4.460E-05, 0, 0, 0, 0, 0, 0, 0, 0], 99645816.200, 3445.050, 516.850, 0.400)
comp0119 = Species(119, 'C13H14', 'Trimethylnaphtalene', 170, 'Mayer-Kelly', [5.36000, 0.15800, -4.980E-05, 0, 0, 0, 0, 0, 0, 0, 0], 66151424.200, 3343.725, 536.850, 0.415)
comp0120 = Species(120, 'C15H18', 'Dimethylpropylnaphthalene', 198, 'Mayer-Kelly', [6.57000, 0.19400, -6.190E-05, 0, 0, 0, 0, 0, 0, 0, 0], 18840595.500, 3039.750, 566.850, 0.410)
comp0121 = Species(121, 'C15H16', 'ДMAH (?)', 196, 'Mayer-Kelly', [8.13000, 0.17900, -5.630E-05, 0, 0, 0, 0, 0, 0, 0, 0], 144444565.500, 3039.750, 566.850, 0.375)
comp0122 = Species(122, 'C10H21*', '', 141, 'Mayer-Kelly', [8.18000, 0.19400, -7.050E-05, 0, 0, 0, 0, 0, 0, 0, 0], -249826296.330, 2107.560, 344.450, 0.490)
comp0123 = Species(123, 'n-C10H19*', '', 139, 'Mayer-Kelly', [9.74000, 0.17900, -6.500E-05, 0, 0, 0, 0, 0, 0, 0, 0], -124222326.330, 2208.885, 341.850, 0.491)
comp0124 = Species(124, 'i-C10H21*', '', 141, 'Mayer-Kelly', [8.18000, 0.19400, -7.050E-05, 0, 0, 0, 0, 0, 0, 0, 0], -249826296.330, 2107.560, 344.450, 0.490)
comp0125 = Species(125, 'C8H15*', '', 111, 'Mayer-Kelly', [-0.97900, 0.17290, -9.640E-05, 2.072E-08, 0, 0, 0, 0, 0, 0, 0], -82982356.180, 2624.318, 293.450, 0.386)
comp0126 = Species(126, 'C12H15*', '', 159, 'Mayer-Kelly', [-10.34000, 0.22300, -8.890E-05, 0, 0, 0, 0, 0, 0, 0, 0], 95877697.100, 3039.750, 516.850, 0.590)
comp0127 = Species(127, 'C13H13*', '', 169, 'Mayer-Kelly', [5.36000, 0.15800, -4.980E-05, 0, 0, 0, 0, 0, 0, 0, 0], 66151424.200, 3343.725, 536.850, 0.415)
comp0128 = Species(128, 'n-C13H27*', '', 183, 'Mayer-Kelly', [11.51000, 0.24900, -9.110E-05, 0, 0, 0, 0, 0, 0, 0, 0], -320708803.400, 1722.525, 384.650, 0.623)
comp0129 = Species(129, '5-C14H29*', '', 197, 'Mayer-Kelly', [13.03000, 0.26900, -9.970E-05, 0, 0, 0, 0, 0, 0, 0, 0], -349179036.600, 1621.200, 420.850, 0.630)
comp0130 = Species(130, 'n-C17H35*', '', 239, 'Mayer-Kelly', [15.95000, 0.32300, -1.190E-04, 0, 0, 0, 0, 0, 0, 0, 0], -407794222.600, 1317.225, 459.850, 0.770)
comp0131 = Species(131, 'n-C14H27*', '', 195, 'Mayer-Kelly', [14.18000, 0.25300, -9.240E-05, 0, 0, 0, 0, 0, 0, 0, 0], -206660398.640, 1560.405, 415.850, 0.644)
comp0132 = Species(132, 'i-C10H19*', '', 139, 'Mayer-Kelly', [10.15000, 0.18100, -6.670E-05, 0, 0, 0, 0, 0, 0, 0, 0], -136489647.400, 2127.825, 339.850, 0.533)
comp0133 = Species(133, 'i-C17H35*', '', 239, 'Mayer-Kelly', [16.36000, 0.32500, -1.200E-04, 0, 0, 0, 0, 0, 0, 0, 0], -414493101.000, 1317.225, 459.850, 0.760)
comp0134 = Species(134, 'i-C7H13*', '', 97, 'Mayer-Kelly', [-0.78900, 0.15040, -8.390E-05, 1.817E-08, 0, 0, 0, 0, 0, 0, 0], -62341437.110, 2837.100, 264.050, 0.358)
comp0135 = Species(135, 'C11H21*', '', 153, 'Mayer-Kelly', [-2.00500, 0.25170, -1.390E-04, 2.954E-08, 0, 0, 0, 0, 0, 0, 0], -270467215.400, 1965.705, 365.650, 0.535)
comp0136 = Species(136, 'cyC5H7*', '', 67, 'Mayer-Kelly', [3.02000, 0.06780, -1.480E-05, 0, 0, 0, 0, 0, 0, 0, 0], 39774590.500, 4053.000, 226.850, 0.193)
comp0137 = Species(137, 'dicyC9H8', '', 116, 'Mayer-Kelly', [3.10000, 0.10900, -3.420E-05, 0, 0, 0, 0, 0, 0, 0, 0], 162447801.200, 3455.183, 416.850, 0.265)
comp0138 = Species(138, '2-C7H14', '', 98, 'Mayer-Kelly', [-0.78900, 0.15000, -8.390E-05, 0, 0, 0, 0, 0, 0, 0, 0], -62341437.110, 2837.100, 264.050, 0.358)
comp0139 = Species(139, 'C11H11*', '', 143, 'Mayer-Kelly', [3.41000, 0.12700, -4.070E-05, 0, 0, 0, 0, 0, 0, 0, 0], 133140208.200, 3566.640, 498.850, 0.345)
comp0140 = Species(140, 'C9H9*', '', 117, 'Mayer-Kelly', [-7.45000, 0.10400, -3.150E-05, 0, 0, 0, 0, 0, 0, 0, 0], 116393012.200, 3404.520, 390.850, 0.244)
comp0141 = Species(141, 'C12H13*', '', 157, 'Mayer-Kelly', [-10.34000, 0.22300, -8.890E-05, 0, 0, 0, 0, 0, 0, 0, 0], 95877697.100, 3039.750, 516.850, 0.590)
comp0142 = Species(142, 'C15H17*', '', 197, 'Mayer-Kelly', [6.57000, 0.19400, -6.190E-05, 0, 0, 0, 0, 0, 0, 0, 0], 18840595.500, 3039.750, 566.850, 0.410)
comp0143 = Species(143, 'dicyC9H10', '', 118, 'Mayer-Kelly', [3.10000, 0.10900, -3.420E-05, 0, 0, 0, 0, 0, 0, 0, 0], 56521786.500, 3232.268, 407.850, 0.265)
comp0144 = Species(144, 'C15H21*', '', 201, 'Mayer-Kelly', [-10.34000, 0.29500, -1.170E-04, 0, 0, 0, 0, 0, 0, 0, 0], -119868055.370, 3039.750, 566.850, 0.380)
comp0145 = Species(145, 'n-C7H12', '', 96, 'Mayer-Kelly', [0.77100, 0.13600, -7.310E-05, 0, 0, 0, 0, 0, 0, 0, 0], 63262532.890, 4053.000, 256.850, 0.680)
comp0146 = Species(146, 'C15H29*', '', 209, 'Mayer-Kelly', [12.58000, 0.27500, -9.970E-05, 0, 0, 0, 0, 0, 0, 0, 0], -333478540.350, 1459.080, 426.850, 0.660)
comp0147 = Species(147, 'C11H10', 'Methylnaphthalene', 142, 'Mayer-Kelly', [3.41000, 0.12700, -4.070E-05, 0, 0, 0, 0, 0, 0, 0, 0], 133140208.200, 3566.640, 498.850, 0.345)
comp0148 = Species(148, 'C14H14', 'Dimethylvinylnaphthalene', 182, 'Mayer-Kelly', [6.93000, 0.16000, -4.950E-05, 0, 0, 0, 0, 0, 0, 0, 0], 161359233.460, 3039.750, 551.850, 0.491)
comp0149 = Species(149, 'C6H5*', '', 77, 'Mayer-Kelly', [-14.60000, 0.14000, -1.080E-04, 3.300E-08, 0, 0, 0, 0, 0, 0, 0], 82479940.300, 4924.395, 288.850, 0.212)
comp0150 = Species(150, '1-C7H13*', '', 97, 'Mayer-Kelly', [-0.78900, 0.15040, -8.390E-05, 1.817E-08, 0, 0, 0, 0, 0, 0, 0], -62341437.110, 2837.100, 264.050, 0.358)
comp0151 = Species(151, 'n-C8H15*', '', 111, 'Mayer-Kelly', [-0.97900, 0.17290, -9.640E-05, 2.072E-08, 0, 0, 0, 0, 0, 0, 0], -82982356.180, 2624.318, 293.450, 0.386)
comp0152 = Species(152, 'n-C9H20', '', 128, 'Mayer-Kelly', [0.75100, 0.16180, -4.610E-05, -7.120E-09, 0, 0, 0, 0, 0, 0, 0], -229185377.260, 2310.210, 321.450, 0.444)
comp0153 = Species(153, 'n-C9H19*', '', 127, 'Mayer-Kelly', [0.75100, 0.16180, -4.610E-05, -7.120E-09, 0, 0, 0, 0, 0, 0, 0], -229185377.260, 2310.210, 321.450, 0.444)
comp0154 = Species(154, '1-C8H16', '', 112, 'Mayer-Kelly', [-0.97900, 0.17290, -9.640E-05, 2.072E-08, 0, 0, 0, 0, 0, 0, 0], -82982356.180, 2624.318, 293.450, 0.386)
comp0155 = Species(155, 'n-C11H24', '', 156, 'Mayer-Kelly', [-2.00500, 0.25170, -1.390E-04, 2.954E-08, 0, 0, 0, 0, 0, 0, 0], -270467215.400, 1965.705, 365.650, 0.535)
comp0156 = Species(156, 'n-C11H23*', '', 155, 'Mayer-Kelly', [-2.00500, 0.25170, -1.390E-04, 2.954E-08, 0, 0, 0, 0, 0, 0, 0], -270467215.400, 1965.705, 365.650, 0.535)
comp0157 = Species(157, '1-C9H18', '', 126, 'Mayer-Kelly', [-0.88800, 0.19400, -1.070E-04, 2.318E-08, 0, 0, 0, 0, 0, 0, 0], -103581407.260, 2340.608, 318.850, 0.430)
comp0158 = Species(158, 'n-C8H18', '', 114, 'Mayer-Kelly', [-1.45600, 0.18420, -1.000E-04, 2.115E-08, 0, 0, 0, 0, 0, 0, 0], -208586326.180, 2482.463, 295.650, 80.394)
comp0159 = Species(159, 'n-C8H17*', '', 113, 'Mayer-Kelly', [-1.45600, 0.18420, -1.000E-04, 2.115E-08, 0, 0, 0, 0, 0, 0, 0], -208586326.180, 2482.463, 295.650, 80.394)
comp0160 = Species(160, '1-C7H14', '', 98, 'Mayer-Kelly', [-0.78900, 0.15040, -8.390E-05, 1.817E-08, 0, 0, 0, 0, 0, 0, 0], -62341437.110, 2837.100, 264.050, 0.358)
comp0161 = Species(161, '1-C9H17', '', 125, 'Mayer-Kelly', [-0.88800, 0.19400, -1.070E-04, 2.318E-08, 0, 0, 0, 0, 0, 0, 0], -103581407.260, 2340.608, 318.850, 0.430)
comp0162 = Species(162, 'O2', '', 32, 'Mayer-Kelly', [6.71300, 0.00000, 4.170E-06, -2.540E-09, 0, 0, 0, 0, 0, 0, 0], 0.000, 5045.985, -118.550, 0.021)
comp0163 = Species(163, 'C2H2O', '', 42, 'Mayer-Kelly', [1.52530, 0.00913, -2.590E-05, 6.445E-09, 0, 0, 0, 0, 0, 0, 0], -61127265.400, 6484.800, 106.850, 0.207)
comp0164 = Species(164, 'CH3COH', '', 44, 'Mayer-Kelly', [1.84340, 0.00353, -2.400E-05, 5.685E-09, 0, 0, 0, 0, 0, 0, 0], -166467128.240, 5572.875, 187.850, 0.303)
comp0165 = Species(165, 'HO2*', '', 33, 'Mayer-Kelly', [4.25000, 0.02690, -2.500E-05, 8.850E-09, 0, 0, 0, 0, 0, 0, 0], -136070967.500, 15198.750, 326.850, 0.090)
comp0166 = Species(166, 'CH2**', '', 14, 'Mayer-Kelly', [0, 0.02520, -9.200E-06, 1.170E-09, 0, 0, 0, 0, 0, 0, 0], -74525022.200, 4640.685, -82.150, 0.008)
comp0167 = Species(167, 'CH3COH*', '', 44, 'Mayer-Kelly', [1.84340, 0.00353, -2.400E-05, 5.685E-09, 0, 0, 0, 0, 0, 0, 0], -166467128.240, 5572.875, 187.850, 0.303)
comp0168 = Species(168, 'CHCO*', '', 41, 'Mayer-Kelly', [1.52530, 0.00913, -2.590E-05, 6.445E-09, 0, 0, 0, 0, 0, 0, 0], -61127265.400, 6484.800, 106.850, 0.207)
comp0169 = Species(169, 'CH***', '', 13, 'Mayer-Kelly', [0, 0.02520, -9.200E-06, 1.170E-09, 0, 0, 0, 0, 0, 0, 0], -74525022.200, 4640.685, -82.150, 0.008)
comp0170 = Species(170, 'C2H*', '', 25, 'Mayer-Kelly', [-2.55000, 0.04950, -3.500E-05, 1.040E-08, 0, 0, 0, 0, 0, 0, 0], 52334987.500, 5066.250, 8.850, 0.065)
comp0171 = Species(171, 'N2', '', 28, 'Mayer-Kelly', [7.44000, -0.00324, 6.400E-06, -2.540E-09, 0, 0, 0, 0, 0, 0, 0], 0.000, 3394.388, -146.950, 0.040)


compset1 = [comp0001,  comp0002,  comp0003,  comp0004,  comp0005,  comp0006,  comp0007,  comp0008,  comp0009,  comp0010,  comp0011,  comp0012,  comp0013,  comp0014,  comp0015,  comp0016,  comp0017,  comp0018,  comp0019,  comp0020,  comp0021,  comp0022,  comp0023,  comp0024,  comp0025,  comp0026,  comp0027,  comp0028,  comp0029,  comp0030,  comp0031,  comp0032,  comp0033,  comp0034,  comp0035,  comp0036,  comp0037,  comp0038,  comp0039,  comp0040,  comp0041,  comp0042,  comp0043,  comp0044,  comp0045,  comp0046,  comp0047,  comp0048,  comp0049,  comp0050,  comp0051,  comp0052,  comp0053,  comp0054,  comp0055,  comp0056,  comp0057,  comp0058,  comp0059,  comp0060,  comp0061,  comp0062,  comp0063,  comp0064,  comp0065,  comp0066,  comp0067,  comp0068,  comp0069,  comp0070,  comp0071,  comp0072,  comp0073,  comp0074,  comp0075,  comp0076,  comp0077,  comp0078,  comp0079,  comp0080,  comp0081,  comp0082,  comp0083,  comp0084,  comp0085,  comp0086,  comp0087,  comp0088,  comp0089,  comp0090,  comp0091,  comp0092,  comp0093,  comp0094,  comp0095,  comp0096,  comp0097,  comp0098,  comp0099,  comp0100,  comp0101,  comp0102,  comp0103,  comp0104,  comp0105,  comp0106,  comp0107,  comp0108,  comp0109,  comp0110,  comp0111,  comp0112,  comp0113,  comp0114,  comp0115,  comp0116,  comp0117,  comp0118,  comp0119,  comp0120,  comp0121,  comp0122,  comp0123,  comp0124,  comp0125,  comp0126,  comp0127,  comp0128,  comp0129,  comp0130,  comp0131,  comp0132,  comp0133,  comp0134,  comp0135,  comp0136,  comp0137,  comp0138,  comp0139,  comp0140,  comp0141,  comp0142,  comp0143,  comp0144,  comp0145,  comp0146,  comp0147,  comp0148,  comp0149,  comp0150,  comp0151,  comp0152,  comp0153,  comp0154,  comp0155,  comp0156,  comp0157,  comp0158,  comp0159,  comp0160,  comp0161,  comp0162,  comp0163,  comp0164,  comp0165,  comp0166,  comp0167,  comp0168,  comp0169,  comp0170,  comp0171]