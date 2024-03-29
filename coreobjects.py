"""
Header      : core_objects.py
Created     : 08.01.2023
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Contains core objects which are base for all modules

References	:

"""

import sys
import numpy as np
import databases.PRBKV_database as binary_db

class UnitsConverter:
	"""
		Describes units conversions for pressure, flowrate and temperature units

		Classes
		----------
		- .Pressure(pressure: float)
			Converts pressure units
		- .Flowrate(flowrate: float)
			Converts flowrate units
		- .Temperature(temperature: float)
			Converts temperature units
		"""
	class Pressure:
		def bara_to_kgfpcm2g(pressure_bara: float):
			return pressure_bara * 1.0197162129779282 - 1.033227
		def kgfpcm2g_to_bara(pressure_kgfpcm2: float):
			return (pressure_kgfpcm2 + 1.033227) / 1.0197162129779282
		def psi_to_kgfpcm2(pressure_psi: float):
			return pressure_psi * 0.0703069579640175
		def bar_to_psi(pressure_bar: float):
			return pressure_bar * 14.503773773
		def bar_to_kPa(pressure_bar: float):
			return pressure_bar * 100
		def kPa_to_psi(pressure_kPa: float):
			return pressure_kPa * 0.1450377377
		def MPa_to_psi(pressure_MPa: float):
			return pressure_MPa * 1000 * 0.1450377377
		def psi_to_kPa(pressure_psi: float):
			return pressure_psi / 0.1450377377
	class Flowrate:
		def sm3d_to_sm3y(flowrate_sm3pday: float):
			return flowrate_sm3pday * 365
	class Temperature:
		def C_to_R(temperature_C: float):
			return (temperature_C + 273.15) * 9 / 5
		def R_to_K(temperature_R: float):
			return temperature_R * 5 / 9
		def C_to_K(temperature_C: float):
			return temperature_C + 273.15
		def K_to_R(temperature_K: float):
			return temperature_K * 1.8


class Species:
    '''
    Describes chemical species

    Methods
    ----------
    .CPIG(T: float)
        Calculates Heat Capacity (constant Pressure) at specified temperature
    .HIGV(T: float)
        Calculates Pure Component Ideal Gas Enthalpy of vapor phase
    '''
    def __init__(self, ID: int, formula: str, name: str, MW: float, CPIGoption: str, CPIGcoeffs: list[float], DHFORM: float,
                 PC: float, TC: float, OMEGA: float):
        '''
        :param ID: Species ID
        :param name: Species Name
        :param: formula: Species Formula
        :param MW: [g/mol] Molar Weight
        :param CPIGoption: Method for CPIG calculation
            1 - DIPPR Equation 107 Heat Capacity Correlation
            2 - Mayer-Kelly Heat Capacity Equation
        :param CPIGcoeffs: Coefficients for Heat Capacity equation (list of 7 coeffs)
        :param DHFORM: [J/kgmol] Pure Component Ideal Gas Enthalpy of Formation @ 25 degC
        :param PC: [kPa] Critical Pressure
        :param TC: [C] Critical Pressure
        :param OMEGA: [dmls.] Pitzer Acentric Factor
        '''
        self.ID = ID
        self.name = name
        self.formula = formula
        self.MW = MW
        self.CPIGoption = CPIGoption
        self.CPIGcoeffs = CPIGcoeffs
        self.DHFORM = DHFORM
        self.PC = PC
        self.TC = TC
        self.OMEGA = OMEGA

    def CPIG(self, T: float):
        '''
        Returns Specific Heat Capacity [J/(mol*K)] of pure component at specified Temperature.
        Two options for calculation are available and must be specified in CPIGoption variable:
            - 'DIPPR eq 107' - DIPPR Equation 107 Heat Capacity Correlation
                coefficients list format: [C1, C2, C3, C4, C5, C6, C7];
            - 'Mayer-Kelly' - Mayer-Kelly Heat Capacity Equation
                coefficients list format: [a, b, c, d, 0, 0, 0]
                WARNING: from available database applicability limits are not known;
            - 'Polynomial' - Aspen Property Ideal Gas heat capacity polynomial
                coefficients format: [C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11]

        :param T: [K] Temperature
        :return: [J/(mol*K)] Specific Heat Capacity
        '''
        CPIG = -1  # Initial value (negative to track malfunctions)
        if self.CPIGoption == 'DIPPR eq 107':
            C1 = self.CPIGcoeffs[0]
            C2 = self.CPIGcoeffs[1]
            C3 = self.CPIGcoeffs[2]
            C4 = self.CPIGcoeffs[3]
            C5 = self.CPIGcoeffs[4]
            C6 = self.CPIGcoeffs[5]
            C7 = self.CPIGcoeffs[6]
            if C6 <= T <= C7:
                CPIG = 4.1868 * (C1 + C2 * (C3 / T / np.sinh(C3 / T)) ** 2 + C4 * (C5 / T / np.cosh(C5 / T)) ** 2)
            else:
                print('ERROR! DIPPR Equation 107 heat capacity correlation is not suitable for specified temperature'
                      'for component {}'.format(self.name))
                sys.exit()
                #CPIG = 4.1868 * (C1 + C2 * (C3 / T / np.sinh(C3 / T)) ** 2 + C4 * (C5 / T / np.cosh(C5 / T)) ** 2)
        elif self.CPIGoption == 'Mayer-Kelly':
            a = self.CPIGcoeffs[0]
            b = self.CPIGcoeffs[1]
            c = self.CPIGcoeffs[2]
            d = self.CPIGcoeffs[3]
            CPIG = 4.1868 * (a + b * T + c / T**2 + d * T**2)

        elif self.CPIGoption == 'Polynomial':
            C1 = self.CPIGcoeffs[0]
            C2 = self.CPIGcoeffs[1]
            C3 = self.CPIGcoeffs[2]
            C4 = self.CPIGcoeffs[3]
            C5 = self.CPIGcoeffs[4]
            C6 = self.CPIGcoeffs[5]
            C7 = self.CPIGcoeffs[6]
            C8 = self.CPIGcoeffs[7]
            C9 = self.CPIGcoeffs[8]
            C10 = self.CPIGcoeffs[9]
            C11 = self.CPIGcoeffs[10]
            if C7 <= T <= C8:
                CPIG = 4.1868 * (C1 + C2 * T + C3 * T **2 + C4 * T**3 + C5 * T**4 + C6 * T**5)
            elif T < C7:
                CPIG = C9 + C10 * T**C11
            else:
                print('ERROR! Aspen Properties polynomial heat capacity correlation is not suitable for specified'
                      'temperature for component {}'.format(self.name))
                sys.exit()
        else:
            print('ERROR! Selected Heat Capacity calculation method for component {} is not available. Specify valid'
              'heat capacity method or check spelling'.format(self.name))
        return CPIG

    def HIGV(self, T: float, numsteps= 10000):
        '''
        Returns Pure Component Ideal Gas Enthalpy of vapor phase [J/kgmol] calculated by direct integration
        (Results for propane converge with Aspen Plus to forth digit)

        :param T: [K] Temperature
        :return: [J/kgmol] Pure Component Ideal Gas Enthalpy
        '''
        init_T = 25 + 273.15
        dT = (T - init_T) / numsteps
        T_vector = list(map(lambda x: init_T + dT * x, range(numsteps + 1)))
        dH_vector = np.array(list(map(lambda x: dT * self.CPIG(x), T_vector)))
        H = self.DHFORM + dH_vector.sum() * 1000
        # Reduce computing time!
        # Use Barin Equation for Entahlpy (in HYSYS: Component --> TDep)
        return H


class Stream:
    '''
        Describes material stream (IDEAL property method?)

        Attributes
        ----------
        - MW [kg/kgmol] Molar weight
        - RHOIG [kg/m3] Ideal gas density @ Actual conditions
        - STDRHOIG [kg/sm3] Ideal gas density @ Std. Conditions
        - CPIG [kJ/(kgmol*K)] Stream ideal gas heat capacity at constant pressure
        - FLMASS [kg/hr] Mass flow
        - FLVOLIG [m3/hr] Actual ideal gas volume flow
        - FLSTDVOLIG [sm3/h3]  Standard ideal gas volume flow
        - COMPMASSFR [mass. frac.] Stream composition in terms of mass fractions
        - COMPMOLCIG [kgmol/m3] Stream composition in terms of molar concentrations (ideal gas)
        '''
    def __init__(self, compset: list[Species], x: dict, molflow: float, P: float, T: float, eos_option: str):
        '''
        :param compset: Set of components in stream
        :param x: [mol. fract.] Stream composition in terms of molar fractions
        :param molflow: [kgmol/hr] Stream molar flow
        :param P: [MPa] Stream pressure
        :param T: [K] Steam temperature
        :param eos_option: Allows to select desired equation of state for VLE calculations. Available options:
            - 'IG'
            - 'PENG-ROB'
        '''
        self.compset = compset
        comp_keys = list(map(lambda x: x.ID, self.compset))  # Keys for component-dependent attributes dictionaries
        if abs(sum(x.values()) - 1 < 1e-2):  # Check for sum of all comps fractions
            if len(x) < len(comp_keys):  # Allows to enter only non-zero comps to Stream input
                x_dict = dict()
                for comp in comp_keys:
                    if comp in x.keys():
                        x_dict[comp] = x[comp]
                    else:
                        x_dict[comp] = 0
                self.COMPMOLFR = x_dict
            else:
                self.COMPMOLFR = x
        else:
            print('\nERROR! Components concentrations for stream #PLACEYORSTREAMNAMEHERE# are not entered correctly '
                  '(sum <>1). Check concentrations input')
            sys.exit()
        self.FLMOL = molflow
        self.P = P
        self.T = T
        self.eos_option = eos_option
        R = 8.31446261815324  # [J/(mole*K)] Gas Constant
        R_field = 10.731577089016  # [psi*ft3/(lbmol*R)] Gas Constant

        '''[kg/kgmol] Stream molar weight'''
        self.MW = sum(list(map(lambda x: self.COMPMOLFR[x.ID] * x.MW, self.compset)))

        '''[kgmol/hr] Molar flow of individual components'''
        indmolflows = list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.ID], self.compset))
        self.FLINDMOL = dict(zip(comp_keys, indmolflows))

        '''[mass. fract.] Stream composition in terms of mass fractions'''
        x_mass = list(map(lambda x: self.COMPMOLFR[x.ID] * x.MW / self.MW, self.compset))
        self.COMPMASSFR = dict(zip(comp_keys, x_mass))

        '''[kg/hr] Mass flow of individual components'''
        indmassflows = list(map(lambda x: self.FLINDMOL[x.ID] * x.MW, self.compset))
        self.FLINDMASS = dict(zip(comp_keys, indmassflows))

        '''[J/(mol*K)] Stream ideal gas heat capacity at constant pressure'''
        self.CP = sum(list(map(lambda x: self.COMPMOLFR[x.ID] * x.CPIG(T), self.compset)))
        # for comp in self.compset:
        #     print('COMPMOLFR', self.COMPMOLFR[comp.name])
        #     print('CPIG', comp.CPIG(T))
        #     print('ID', comp.ID)
        #     print('T', T)
        #     print()
        # Check convergence with Hysys! - It does not converge with Hysys.

        if self.eos_option == 'IG':
            '''[kg/m3] Stream density @ Actual conditions (Ideal Gas)'''
            self.RHO = self.MW / (R * self.T / (self.P * 1000))

            '''[kg/sm3] Stream density @ Standard conditions'''
            self.STDRHO = self.MW/ (R * self.T)

            '''[kg/hr] Stream mass flow'''
            self.FLMASS = sum(list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.ID] * x.MW, self.compset)))

            '''[m3/hr] Stream volume flow @ Actual conditions (Ideal Gas)'''
            self.FLVOL = self.FLMOL * self.MW / self.RHO

            '''[sm3/hr] Stream volume flow @ Standard conditions'''
            self.STDVOL = self.FLMOL * self.MW/ self.STDRHO

            '''[kgmol/m3] Stream composition in terms of molar concentrations @ Actual Conditions (Ideal Gas)'''
            molconc = list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.ID] / self.FLVOL, self.compset))
            self.COMPMOLC = dict(zip(comp_keys, molconc))

        elif self.eos_option == 'PENG-ROB':
            PRKBV1_df = binary_db.PRKBV1.loc[comp_keys, comp_keys]
            '''Make sure that order of elements is same as for other arrays'''
            PRKBV1_df = PRKBV1_df.reindex(columns=comp_keys, index=comp_keys)
            # WARNING! If there are components in keys list that are not present in database DataFrame, NaN cells would be created
            '''Z-factor'''
            # FOR ONE PHASE (VAPOR) ONLY!
            '''Component-dependent variables'''
            Pc_arr = np.array(list(map(lambda x: UnitsConverter.Pressure.kPa_to_psi(x.PC), self.compset)))  # Array of PC [psi]
            Tc_arr = np.array(list(map(lambda x: UnitsConverter.Temperature.C_to_R(x.TC), self.compset)))  # Array of TC [R]
            w_arr = np.array(list((map(lambda x: x.OMEGA, self.compset)))) # Array of Acentric factors [dmls.]
            Tr_arr = UnitsConverter.Temperature.K_to_R(self.T) / Tc_arr  # Array of Reduced Temperatures [R/R]
            kappa_arr = np.where(w_arr > 0.49,
                                 0.379642 + 1.4853 * w_arr - 0.164423 * w_arr ** 2 + 0.01666 * w_arr ** 3,
                                 0.37464 + 1.5422 * w_arr - 0.26992 * (w_arr ** 2))
            '''For heavy hydrocarbon comps (w>0.49) used 1980's modification'''
            alfa_arr = np.where(np.logical_and(Tc_arr == 374.149011230469, Tr_arr ** 0.5 < 0.85),  # why second condition?
                            (1.0085677 + 0.82154 * (1 - Tr_arr ** 0.5)) ** 2,
                            (1 + kappa_arr * (1 - Tr_arr ** 0.5)) ** 2)
            '''For water used 1980's modification'''
            ac_arr = 0.45724 * (R_field ** 2) * (Tc_arr ** 2) / Pc_arr
            ai_arr = ac_arr * alfa_arr  # Array of first comp-dependent variables 'ai'
            bi_arr = 0.07780 * R_field * Tc_arr / Pc_arr  # Array of second comp-dependent variables 'bi'

            ''' Mixing rules (eq. 3-22 from [1])'''
            xij_arr = np.array(list(map(lambda x: self.COMPMOLFR[x], comp_keys)))  # Array of comp molar fractions [mol. fract.]
            xlj_arr = np.reshape(xij_arr, (len(xij_arr), 1))  # Array of comp molar fractions [mol. fract.] in columns
            al_arr = np.reshape(ai_arr, (len(ai_arr), 1))  # Array of first comp-dependent variables 'ai' in columns
            deltail_matr = np.array(PRKBV1_df)  # Matrix of binary interaction coefficients 'deltail'
            aj = np.sum(xij_arr * xlj_arr * (ai_arr * al_arr) ** 0.5 * (1 - deltail_matr))  # First mixing rule variable 'aj'
            bj = np.sum(xij_arr * bi_arr)  # Second mixing rule variable 'bj'

            '''Phase-dependent variables (eq. 3-29, 3-30 from [1])'''  # FOR ONE PHASE (VAPOR) ONLY!
            P_field = UnitsConverter.Pressure.MPa_to_psi(self.P)
            T_field = UnitsConverter.Temperature.K_to_R(self.T)
            Aj = aj * P_field / R_field ** 2 / T_field ** 2  # First phase-dependent variable 'Aj'
            Bj = bj * P_field / R_field / T_field  # Second phase-dependent variable 'Bj'

            '''Phase-component-dependent variables for fugacities calculations (eq. 3-20 an 3-21 from [1])'''
            # FOR ONE PHASE (VAPOR) ONLY!
            Aijprime_arr = 1 / aj * 2 * ai_arr ** 0.5 * np.sum(xij_arr * ai_arr ** 0.5 * (1 - deltail_matr), axis= 1)  # Array of first phase-component-dependent variables 'Aij''
            Bijprime = bi_arr / bj

            '''Z-factor (eq. 3-31 from [1], solution of cubic equation from [2])'''
            a2 = -(1 - Bj)
            a1 = Aj - 2 * Bj - 3 * Bj**2
            a0 = -(Aj * Bj - Bj ** 2 - Bj ** 3)
            q = 1 / 3 * a1 - 1 / 9 * a2 ** 2
            r = 1 / 6 * (a1 * a2 - 3 * a0) - 1 / 27 * a2 ** 3
            ''' Description of possible roots:
            if (q**3 + r**2) > 0:
                Equation has one real root and a pair of complex conjugate roots...
            elif (q**3 + r**2) == 0:
                Equation has all real roots and at leas two of them are equal...
            else:
                Equation has all real roots...'''
            '''smallvar introduced to avoid problems when expression is too small '''
            if abs(q ** 3 + r ** 2) <= 1e-8:  # this convergence criteria may affect equicomp results
                smallvar = 0
            else:
                smallvar = (q ** 3 + r ** 2) ** (1 / 2)
            s1 = np.cbrt(r + smallvar)
            s2 = np.cbrt(r - smallvar)
            z1 = (s1 + s2) - a2 / 3
            z2 = complex(-1 / 2 * (s1 + s2) - a2 / 3, (3 ** (1 / 2)) / 2 * (s1 - s2))
            z3 = complex(-1 / 2 * (s1 + s2) - a2 / 3, -(3 ** (1 / 2)) / 2 * (s1 - s2))
            '''Check roots as in 3.7.33 [2]'''
            check1 = abs((z1 + z2 + z3) - complex(-1 * a2, 0)) < 0.001
            check2 = abs((z1 * z2 + z1 * z3 + z2 * z3) - complex(a1, 0)) < 0.001
            check3 = abs((z1 * z2 * z3) - complex(-1 * a0, 0)) < 0.001
            if check1 and check2 and check3:
                pass
            else:
                print('\nWARNING: Roots for z-factor equations are NOT CHECKED successfuly!\n(try to relax "smallvar" critera)')
            roots = list()
            '''Eliminate small imaginary parts (option to improve stability)'''
            for root in [z1, z2, z3]:
                if abs(root.imag) < 10 ** -6:
                    root = root.real
                roots.append(root)
            '''Choose appropriate for specific phase z-factor (currently only vapor phase)
               Z-factor for vapor phase is the the real one with the highest value'''
            zfactor = max([root for root in roots if (not (type(root) == complex) and root >= 0 and root >= Bj)])
            self.Z = zfactor

            '''[kg/m3] Stream density @ Actual conditions (Peng-Robinson EOS)'''
            self.RHO = self.MW / self.Z / (R * self.T / (self.P * 1000))

            '''[m3/hr] Stream volume flow @ Actual conditions (Peng-Robinson EOS)'''
            self.FLVOL = self.FLMOL * self.MW / self.RHO

            '''[kgmol/m3] Stream composition in terms of molar concentrations @ Actual Conditions (Peng-Robinson EOS)'''
            molconcPR = list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.ID] / self.FLVOL, self.compset))
            self.COMPMOLC = dict(zip(comp_keys, molconcPR))
        else:
            print('\nERROR! Selected EOS for #PLACEYORSTRAMNAMEHERE# is not available. Specify valid EOS or check spelling')
            sys.exit()

    def setcomposition(self, comp_x0: dict[str: float]):

        return 1