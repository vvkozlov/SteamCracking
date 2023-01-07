'''
Header      : rctrs_engine.py
Created     : 09.07.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Mathematical solver for Plug-Flow Reactor

References:
        [1] A. Jebarjadi - Multi-Phase Multi-Component Equilibrium Flash Calculations for CompFlow Bio using
            Modified Volume-Translated Peng-Robinson EOS (2017)'
        [2] M. Abramowitz - Handbook of Mathematical Functions with Formulas, Graphs and Matematical Tables

'''


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import sys
import PRBKV_database as binary_db  # move to separate database


'''One step integration with Euler method for differential equations'''
def increment_euler(f, x0, y0, h):
    return y0 + h * f(x0, y0)


'''One step integration with 4th order Runge-Kutta method for differential equations'''
def increment_rungekutta4th(f, x0, y0, h):
    k1 = f(x0, y0)
    k2 = f(x0 + h / 2, y0 + k1 * h / 2)
    k3 = f(x0 + h / 2, y0 + k2 * h / 2)
    k4 = f(x0 + h, y0 + k3 * h)
    return y0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


'''One step integration with 5th order Runge-Kutta-Felberg method for differential equations'''
def increment_rungekuttafelberg5th(f, x0, y0, h):
    k1 = h * f(x0,                  y0)
    k2 = h * f(x0 + 1 / 4 * h,      y0 + 1 / 4 * k1)
    k3 = h * f(x0 + 3 / 8 * h,      y0 + 3 / 32 * k1 + 9 / 32 * k2)
    k4 = h * f(x0 + 12 / 13 * h,    y0 + 1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3)
    k5 = h * f(x0 + h,              y0 + 439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4)
    k6 = h * f(x0 + 1 / 2 * h,      y0 - 8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5)
    return y0 + 16 / 135 * k1 + 6656 / 12825 * k3 + 28561 / 56430 * k4 - 9 / 50 * k5 + 2 / 55 * k6


class UnitsConverter:
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
            #ATTENTION!
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
    def __init__(self, ID: int, name: str, formula: str, MW: float, CPIGDP: list[float], DHFORM: float,
                 PC: float, TC: float, OMEGA: float):
        '''
        :param ID: Species ID
        :param name: Species Name
        :param: formula: Species Formula
        :param MW: [g/mol] Molar Weight
        :param CPIGDP: [K, cal/mol] Coefficients for DIPPR Equation 107
        :param DHFORM: [J/kgmol] Pure Component Ideal Gas Enthalpy of Formation @ 25 degC
        :param PC: [kPa] Critical Pressure
        :param TC: [C] Critical Pressure
        :param OMEGA: [dmls.] Pitzer Acentric Factor
        '''
        self.ID = ID
        self.name = name
        self.formula = formula
        self.MW = MW
        self.CPIGDP = CPIGDP
        self.DHFORM = DHFORM
        self.PC = PC
        self.TC = TC
        self.OMEGA = OMEGA

    def CPIG(self, T: float):
        '''
        Returns Specific Heat Capacity [J/(mol*K)] of pure component at specified Temperature calculated with
        DIPPR Equation 107 Heat Capacity correlation

        :param T: [K] Temperature
        :return: [J/(mol*K)] Specific Heat Capacity

        previously used fourth-order equation:
        4.1887 * (self.Cp_coeffs[0] + self.Cp_coeffs[1] * T
                         + self.Cp_coeffs[2] * T ** 2 + self.Cp_coeffs[3] * T ** 3) '''
        C1 = self.CPIGDP[0]
        C2 = self.CPIGDP[1]
        C3 = self.CPIGDP[2]
        C4 = self.CPIGDP[3]
        C5 = self.CPIGDP[4]
        C6 = self.CPIGDP[5]
        C7 = self.CPIGDP[6]
        if C6 <= T <= C7:
            CPIG = 4.1868 *(C1 + C2 * (C3 / T / np.sinh(C3 / T)) ** 2 + C4 * (C5 / T / np.cosh(C5 / T)) ** 2)
        else:
            print('ERROR! DIPPR Equation 107 heat capacity correlation is not suitable for specified temperature')
            sys.exit()
            #CPIG = 4.1868 * (C1 + C2 * (C3 / T / np.sinh(C3 / T)) ** 2 + C4 * (C5 / T / np.cosh(C5 / T)) ** 2)
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
        self.COMPMOLFR = x
        self.FLMOL = molflow
        self.P = P
        self.T = T
        R = 8.31446261815324  # [J/(mole*K)] Gas Constant
        R_field = 10.731577089016  # [psi*ft3/(lbmol*R)] Gas Constant
        comp_keys = list(map(lambda x: x.name, self.compset))  # Keys for component-dependent attributes dictionaries
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
        self.RHOPR = self.MW / self.Z / (R * self.T / (self.P * 1000))

        '''[m3/hr] Stream volume flow @ Actual conditions (Peng-Robinson EOS)'''
        self.FLVOLPR = self.FLMOL * self.MW / self.RHOPR

        '''[kgmol/m3] Stream composition in terms of molar concentrations @ Actual Conditions (Peng-Robinson EOS)'''
        molconcPR = list(map(lambda x: self.FLMOL * self.COMPMOLFR[x.name] / self.FLVOLPR, self.compset))
        self.COMPMOLCPR = dict(zip(comp_keys, molconcPR))


class Reaction:
    '''
    Describes kinetic reactions

    Methods
    ----------
    .rate(T: float)
        Calculates Reaction Rate at specified temperature with Arrhenius Law
    '''
    def __init__(self, name: str, reagents: list[Species], stoic: list[float], order: list[float], dH: float, k0: float, E0: float):
        '''
        :param name: Reaction Name
        :param reagents: List of all reagents (same order as in equation)
        :param stoic: List of stoichiometric coefficients for listed components (written in same order as reagents).
        Negative values for reagents (left side of equation), positive for products (right side)
        :param order: List of reaction orders by each reactant for listed components (written in same order as reagents).
        Negative values for reagents (left side of equation), positive for products (right side)
        :param dH: [kJ/mol] Heat of Reaction (if equals zero value will be obtained from reagents enthalpy difference)
        :param k0: Reaction Rate Constant (Arrhenius Parameter)
        :param E0: [kJ/mol] Activation Energy
        '''
        self.name = name
        self.reagents = reagents
        self.stoic = dict(zip(list(map(lambda x: x.name, self.reagents)), stoic))
        self.order = dict(zip(list(map(lambda x: x.name, self.reagents)), order))
        self.k0 = k0
        self.E0 = E0

        DHFORM_vect = np.array(list(map(lambda x: x.DHFORM, reagents)))
        if dH == 0:
            self.dH = np.sum(np.array(stoic) * DHFORM_vect) / 10**6
        else:
            self.dH = dH

    def rate(self, T: float, conc: dict):
        '''
        Returns Reaction Rate for forward rxns at specified temperature

        :param T: [K] Temperature
        :param conc: [kmol/m3] Concentrations of components
        :return: Reaction Rate
        '''
        '''
        equation used: v = k * prod([I] ^ i)
            where k = k0 * exp(E0 / R / T)
        '''
        mult = 1  # [kgkol/m3]^n Multiplier that considers contribution of concentrations
        for comp in self.reagents:
            if self.order[comp.name] < 0:
                mult = mult * ((conc[comp.name]) ** abs(self.order[comp.name]))
        # Needs to be revised (use matrix instead of loop)!

        rate = mult * self.k0 * np.exp(-(self.E0 * 1000) / 8.3144 / T)  # [kgmol/(m3*s)]
        '''
        from Dante, 1979 k0 reported in [l / (mol * s)] --> for second-oreder reactions rate is in [kgmol/(m3*s)]
        '''
        return rate


class PFRreactor:
    '''
    Describes adiabatic Plug-Flow Reactor

    Methods
    ----------

    '''
    def __init__(self, length: float, diameter: float, numtubes: float, rxnset: list[Reaction]):
        '''
        :param length: [m] Reactor tube length
        :param diameter: [m] Reactor tube diameter
        :param numtubes: [No] Number of reactor tubes
        :param rxnset: [list of Reaction] Set of reactions occurring in reactor
        '''
        self.length = length
        self.diameter = diameter
        self.numtubes = numtubes
        self.rxnset = rxnset

    # without matrices 1863.49 ms
    def simulation(self, inlet: Stream, dl: float, log: bool) -> tuple[Stream, pd.DataFrame]:
        '''
        Performs  integration along reactor length with Euler method for reactions listed

        :param inlet: [Stream] Reactor inlet stream
        :param step: [m] Integration resolution (step along reactor length)
        :param log: Select if tabular results for each timestep is required
        :return: [Stream] Reactor outlet stream and [pd.DataFrame] Calculations results on each iteration
        '''
        '''Determine conditions at rctr inlet'''
        flow = inlet
        cell_volume = np.pi * (self.diameter ** 2) / 4 * dl  # [m3]
        l = 0  # [m]
        t = 0  # [s]
        '''Df to store calculation results'''
        temp_df = pd.DataFrame()
        '''Keys for components and reactions lists - to make sure that all matrices are uniform'''
        comp_keys = sorted(flow.COMPMOLFR)  # Some more elegant way to create matching list should be found
        rxn_keys = list(map(lambda x: x.name, self.rxnset))
        '''Reaction stoich coefficients matrix [No. rxns x No. comps]'''
        stoic_df = pd.DataFrame(index= rxn_keys, columns= comp_keys)
        '''Reaction order matrix [No. rxns x No. comps]'''
        order_df = pd.DataFrame(index=rxn_keys, columns=comp_keys)
        '''Reaction enthalpy difference dH matrix [No. rxns x 1]'''
        rxndH_df = pd.DataFrame(index=rxn_keys, columns=['dH'])
        '''Assemble stoich coeffs and rxn enthalpies df's'''
        for rxn in self.rxnset:
            rxndH_df['dH'][rxn.name] = rxn.dH
            for comp in inlet.compset:
                if comp in rxn.reagents:
                    stoic_df[comp.name][rxn.name] = rxn.stoic[comp.name]
                else:
                    stoic_df[comp.name][rxn.name] = 0
        '''Convert stoich coeffs and rxn enthalpies df's to matrices'''
        stoic_matrix = np.array(stoic_df)
        rxndH_matrix = np.array(rxndH_df)

        '''Integration through reactor length'''
        while l < self.length:
            '''Residence time for finite rctr cell'''
            dt = cell_volume / flow.FLVOLPR * 3600  # [s]
            print('\tintegration l = {:.3f} m'.format(l + dl))
            print('\t            t = {:.3f} s'.format(t + dt))
            '''Determine conditions at cell inlet'''
            act_C = flow.COMPMOLCPR  # [kgmol/m3]
            act_T = flow.T  # [K]
            act_P = flow.P  # [MPa]
            act_Cp = flow.CPIG  # [kJ/(kg*K)]
            '''Comps concentrations vector [1 x No. comps]'''
            C_vect = np.array(list(dict(sorted(act_C.items())).values()))  # [kgmol/m3]

            '''Functional form of PFReactor mass balance differential equation for integration methods'''
            def concentrations_derivative(x, y, _stoic_matrix= stoic_matrix, _rxnset= self.rxnset, _T= act_T):
                x = 1
                '''Reactions rate constants matrix [No. rxns x 1]'''
                _rateconst_matrix = np.array(list(map(lambda x: x.rate(_T, dict(zip(comp_keys, y))), _rxnset)))  # [kgmol/(m3*s)]
                _rateconst_matrix = np.reshape(_rateconst_matrix, (len(_rateconst_matrix), 1))
                return (_stoic_matrix * _rateconst_matrix).sum(axis= 0)
            '''Comps conc at cell outlet from PFReactor diff equation [1 x No. comps]'''
            C_vect = increment_rungekutta4th(concentrations_derivative, 1, C_vect, dt)  # [kgmol/m3]
            '''Update comps concentration dictionary'''
            act_C = dict(zip(comp_keys, C_vect))

            '''Reactions rate constants matrix [No. rxns x 1] (for new concentrations?)'''
            rateconst_matrix = np.array(list(map(lambda x: x.rate(act_T, act_C), self.rxnset)))  # [kgmol/(m3*s)]
            rateconst_matrix = np.reshape(rateconst_matrix, (len(rateconst_matrix), 1))
            '''Sum of reaction heat for all rxns in rctr [1 x 1]'''
            dQ = np.sum(rateconst_matrix * rxndH_matrix) * -1000  # [kJ/(m3*s)]
            '''Functional form of PFReactor heat balance differential equation for integration methods'''
            def temperature_derivative(x, y, _dQ= dQ, _P= act_P, _Cp= act_Cp):
                x = 1
                return _dQ * 0.008314463 * y / _P / _Cp
            '''Update cell temperature'''
            new_T = increment_rungekutta4th(temperature_derivative, 1, act_T, dt)
            '''Comps mole fractions at cell outlet'''
            new_compmolfr = dict(zip(comp_keys, list(map(lambda x: act_C[x] / sum(act_C.values()), comp_keys))))  # [mol. fract.]
            '''Comps mole flow at cell outlet (volume calculated from PR EOS or IG EOS at cell inlet)'''
            new_molflow = sum(list(map(lambda x: flow.FLVOLPR * act_C[x], comp_keys)))  # [kgmol/hr]
            '''Update flow to cell outlet conditions'''
            flow = Stream(flow.compset, new_compmolfr, new_molflow, act_P, new_T)
            '''Step forward through reactor'''
            l += dl
            t += dt
            '''Dict to store output variables inside the loop before appending to results dataframe'''
            output_line = dict()
            '''Fill output_line'''
            output_line.update(flow.COMPMOLFR.copy())
            output_line['FLMOL [kgmol/hr]'] = flow.FLMOL
            output_line['l [m]'] = l
            output_line['t [s]'] = t
            output_line['T [K]'] = flow.T
            '''Add output line as new index to output df'''
            temp_df = temp_df.append(output_line, ignore_index=True)
        '''Set lengths column as df indexes for result df'''
        temp_df = temp_df.set_index('l [m]', drop=True)

        ''' Possible problems in future:
        1. matrix indexes may be messed up - check comp_keys and rxn_keys when assigning new non-scalar variables
        2.
        '''

        return flow, temp_df


def plot_results(results: pd.DataFrame):
    '''
    Plots results of reactor simulation. Input DataFrame must be in format defined in ReactorModel.Simulation

    :param results: Tabular results of ReactorModel.Simulation
    '''
    fig, axes = plt.subplots(2, 1)
    axes[0].set_title('Composition and Temperature of Reaction Mixture vs. Time')
    axes[0].grid()
    axes[0].set_ylabel('Molar Fraction, [mol. frac.]')
    for param in results.columns[0:-1]:
        if 'rate' not in param and 'FLMASS' not in param and 'FLMOL' not in param and param != 't [s]':
            axes[0].plot(results.index, results[param], label= param)
    axes[0].legend()
    axes[1].set_ylabel('Temperature, [K]')
    axes[1].set_xlabel('Length, [m]')
    axes[1].plot(results.index, results['T [K]'])
    plt.grid()
    plt.show()
    return 1



