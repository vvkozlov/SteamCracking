"""
Header      : chemistry.py (former rctr_engine.py)
Created     : 09.07.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Handles kinetic reactions and solves Plug-Flow Reactor

References:
        [1] A. Jebarjadi - Multi-Phase Multi-Component Equilibrium Flash Calculations for CompFlow Bio using
            Modified Volume-Translated Peng-Robinson EOS (2017)'
        [2] M. Abramowitz - Handbook of Mathematical Functions with Formulas, Graphs and Matematical Tables
        [3] M. Dente - Detailed Prediction of Olefin Yields from Hydrocarbon Pyrolysis Through a Fundamental
            Simulation Model (Spyro), 1979
        [4] Terrasug-2000 User Manual
        [5] Aspen HYSYS User Manual

"""


import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from coreobjects import Species, Stream
import usermath as m


class Reaction:
    """
    Describes kinetic reactions

    Methods
    ----------
    .rate(T: float)
        Calculates Reaction Rate at specified temperature with Arrhenius Law
    """
    def __init__(self, name: str, reagents: list[Species], stoic: list[float], order: list[float], dH: float, k0: float, E0: float):
        """
        :param name: Reaction Name
        :param reagents: List of all reagents (same order as in equation)
        :param stoic: List of stoichiometric coefficients for listed components (written in same order as reagents).
        Negative values for reagents (left side of equation), positive for products (right side)
        :param order: List of reaction orders by each reactant for listed components (written in same order as reagents).
        Negative values for reagents (left side of equation), positive for products (right side)
        :param dH: [kJ/mol] Heat of Reaction (if equals zero value will be obtained from reagents enthalpy difference)
        :param k0: Reaction Rate Constant (Arrhenius Parameter)
        :param E0: [kJ/mol] Activation Energy
        """
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
        - from Dente, 1979 k0 reported in [l/(mol*s)] --> for second-order reactions rate is in [kgmol/(m3*s)]
        - from Terrasug-2000 k0 reported in [ml/(mol*s)] --> for second-order reactions rate is in 0.001 * [kgmol/(m3*s)]
            so in rxn data input used 10^A * 0.001 thus converting to [l/(mol*s)]
        - if does not work - try use k0 in [l/(mol*s)] same as per Dente, 1979 because modules are to close
        '''
        return rate


class PFReactor:
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
    def simulation(self, inlet: Stream, dl: float, eos_option: str, log: bool) -> tuple[Stream, pd.DataFrame]:
        '''
        Performs  integration along reactor length with Euler method for reactions listed

        :param inlet: [Stream] Reactor inlet stream
        :param step: [m] Integration resolution (step along reactor length)
        :param eos_option: Allows to select desired equation of state for VLE calculations. Available options:
            - 'IG'
            - 'PENG-ROB'
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
            '''Calculate volume flowrate trough cell and determine initial concentrations'''
            if eos_option == "IG":
                volflow = flow.FLVOLIG
                act_C = flow.COMPMOLCIG  # [kgmol/m3]
            elif eos_option == "PENG-ROB":
                volflow = flow.FLVOLPR
                act_C = flow.COMPMOLCPR  # [kgmol/m3]
            else:
                print('ERROR! Selected EOS for is not available. Specify valid EOS or check spelling')
                sys.exit()
            '''Residence time for finite rctr cell'''
            dt = cell_volume / volflow * 3600  # [s]
            print('\tintegration l = {:.3f} m'.format(l + dl))
            print('\t            t = {:.3f} s'.format(t + dt))
            '''Determine conditions at cell inlet'''
            act_T = flow.T  # [K]
            act_P = flow.P  # [MPa]
            act_Cp = flow.CPIG  # [kJ/(kg*K)] - only IG option availible
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
            C_vect = m.integrate('rungekutta4th', concentrations_derivative, 1, C_vect, dt)  # [kgmol/m3]
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
            new_T = m.integrate('rungekutta4th', temperature_derivative, 1, act_T, dt)
            '''Comps mole fractions at cell outlet'''
            new_compmolfr = dict(zip(comp_keys, list(map(lambda x: act_C[x] / sum(act_C.values()), comp_keys))))  # [mol. fract.]
            '''Comps mole flow at cell outlet (volume calculated from PR EOS or IG EOS at cell inlet)'''
            if eos_option == "IG":
                new_molflow = sum(list(map(lambda x: flow.FLVOLIG * act_C[x], comp_keys)))  # [kgmol/hr]
            elif eos_option == "PENG-ROB":
                new_molflow = sum(list(map(lambda x: flow.FLVOLPR * act_C[x], comp_keys)))  # [kgmol/hr]
            else:
                print('ERROR! Selected EOS for is not available. Specify valid EOS or check spelling')
                sys.exit()
            '''Update flow to cell outlet conditions'''
            flow = Stream(flow.compset, new_compmolfr, new_molflow, act_P, new_T)
            '''Step forward through reactor'''
            l += dl
            t += dt
            '''Dict to store output variables inside the loop before appending to results dataframe'''
            output_line = dict()
            '''Fill output_line'''
            if eos_option == "IG":
                output_line.update(flow.COMPMOLIG.copy())
            elif eos_option == "PENG-ROB":
                output_line.update(flow.COMPMOLFR.copy())
            else:
                print('ERROR! Selected EOS for is not available. Specify valid EOS or check spelling')
                sys.exit()
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
