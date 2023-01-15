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
from usermath import UnitsConverter as convert

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
    def __init__(self, ID: int, name: str, formula: str, MW: float, CPIGoption: str, CPIGcoeffs: list[float], DHFORM: float,
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
                print('ERROR! DIPPR Equation 107 heat capacity correlation is not suitable for specified temperature')
                sys.exit()
                #CPIG = 4.1868 * (C1 + C2 * (C3 / T / np.sinh(C3 / T)) ** 2 + C4 * (C5 / T / np.cosh(C5 / T)) ** 2)

        elif self.CPIGoption == 'Mayer-Kelly':
            a = self.CPIGcoeffs[0]
            b = self.CPIGcoeffs[1]
            c = self.CPIGcoeffs[2]
            d = self.CPIGcoeffs[3]
            CPIG = 4.1868 * (a + b * T + c / T**2 + d * T**2)
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
    def __init__(self, compset: list[Species], x: dict, molflow: float, P: float, T: float):
        '''
        :param compset: Set of components in stream
        :param x: [mol. fract.] Stream composition in terms of molar fractions
        :param molflow: [kgmol/hr] Stream molar flow
        :param P: [MPa] Stream pressure
        :param T: [K] Steam temperature
        '''
        self.compset = compset
        comp_keys = list(map(lambda x: x.name, self.compset))  # Keys for component-dependent attributes dictionaries
        if sum(x.values()) == 1:
            if len(x) < len(comp_keys):
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
            print('ERROR! Components concentrations for stream #PLACEYORSTREAMNAMEHERE# are not entered correctly '
                  '(sum <>1). Check concentrations input')
            sys.exit()
        self.FLMOL = molflow
        self.P = P
        self.T = T
        R = 8.31446261815324  # [J/(mole*K)] Gas Constant
        R_field = 10.731577089016  # [psi*ft3/(lbmol*R)] Gas Constant
        PRKBV1_df = binary_db.PRKBV1.loc[comp_keys, comp_keys]
        '''Make sure that order of elements is same as for other arrays'''
        PRKBV1_df = PRKBV1_df.reindex(columns= comp_keys, index= comp_keys)
        # WARNING! If there are components in keys list that are not present in database DataFrame, NaN cells would be created

        '''[kg/kgmol] Stream molar weight'''
        self.MW = sum(list(map(lambda x: self.COMPMOLFR[x.name] * x.MW, self.compset)))

        '''[kg/m3] Stream density @ Actual conditions (Ideal Gas)'''
        self.RHOIG = self.MW / (R * self.T / (self.P * 1000))

        '''[kg/sm3] Stream density @ Standard conditions'''
        self.STDRHOIG = self.MW/ (R * self.T)

        '''[J/(mol*K)] Stream ideal gas heat capacity at constant pressure'''
        self.CPIG = sum(list(map(lambda x: self.COMPMOLFR[x.name] * x.CPIG(T), self.compset)))
        # Check convergence with Hysys! - It does not converge with Hysys.

        '''[kg/hr] Stream mass flow'''
        self.FLMASS = sum(list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.name] * x.MW, self.compset)))

        '''[m3/hr] Stream volume flow @ Actual conditions (Ideal Gas)'''
        self.FLVOLIG = self.FLMOL * self.MW / self.RHOIG

        '''[sm3/hr] Stream volume flow @ Standard conditions'''
        self.STDVOLIG = self.FLMOL * self.MW/ self.STDRHOIG

        '''[kgmol/hr] Molar flow of individual components'''
        indmolflows = list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.name], self.compset))
        self.FLINDMOL = dict(zip(comp_keys, indmolflows))

        '''[mass. fract.] Stream composition in terms of mass fractions'''
        x_mass = list(map(lambda x: self.COMPMOLFR[x.name] * x.MW / self.MW, self.compset))
        self.COMPMASSFR = dict(zip(comp_keys, x_mass))

        '''[kgmol/m3] Stream composition in terms of molar concentrations @ Actual Conditions (Ideal Gas)'''
        molconc = list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.name] / self.FLVOLIG, self.compset))
        self.COMPMOLCIG = dict(zip(comp_keys, molconc))

        '''[kg/hr] Mass flow of individual components'''
        indmassflows = list(map(lambda x: self.FLINDMOL[x.name] * x.MW, self.compset))
        self.FLINDMASS = dict(zip(comp_keys, indmassflows))

        '''Z-factor'''
        # FOR ONE PHASE (VAPOR) ONLY!
        '''Component-dependent variables'''
        Pc_arr = np.array(list(map(lambda x: convert.Pressure.kPa_to_psi(x.PC), self.compset)))  # Array of PC [psi]
        Tc_arr = np.array(list(map(lambda x: convert.Temperature.C_to_R(x.TC), self.compset)))  # Array of TC [R]
        w_arr = np.array(list((map(lambda x: x.OMEGA, self.compset)))) # Array of Acentric factors [dmls.]
        Tr_arr = convert.Temperature.K_to_R(self.T) / Tc_arr  # Array of Reduced Temperatures [R/R]
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
        P_field = convert.Pressure.MPa_to_psi(self.P)
        T_field = convert.Temperature.K_to_R(self.T)
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
        self.RHOPR = self.MW / self.Z / (R * self.T / (self.P * 1000))

        '''[m3/hr] Stream volume flow @ Actual conditions (Peng-Robinson EOS)'''
        self.FLVOLPR = self.FLMOL * self.MW / self.RHOPR

        '''[kgmol/m3] Stream composition in terms of molar concentrations @ Actual Conditions (Peng-Robinson EOS)'''
        molconcPR = list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.name] / self.FLVOLPR, self.compset))
        self.COMPMOLCPR = dict(zip(comp_keys, molconcPR))

    def setcomposition(self, comp_x0: dict[str: float]):

        return 1