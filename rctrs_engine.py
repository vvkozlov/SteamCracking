'''
Header      : rctrs_engine_v1.py
Created     : 09.07.2022
Modified    : 09.07.2022
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Retrieved from reactors_v3.py as individual module
Changes     :
'''
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math


class Species:
    '''
    Describes chemical species

    Methods
    ----------
    .Cp(T: float)
        Calculates Heat Capacity (constant Pressure) at specified temperature
    '''
    def __init__(self, name: str, MW: float, CPIGDP: list, DHFORM: float):
        '''
        :param name: Species Name
        :param MW: [g/mol] Molar Weight
        :param CPIGDP: [K, cal/mol] Coefficients for DIPPR Equation 107
        :param DHFORM: [J/kgmol] Pure Component Ideal Gas Enthalpy of Formation @ 25 degC
        '''
        self.name = name
        self.MW = MW
        self.CPIGDP = CPIGDP
        self.DHFORM = DHFORM

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
            # print('ERROR! DIPPR Equation 107 heat capacity correlation is not suitable for specified temperature')
            # sys.exit()
            CPIG = 4.1868 * (C1 + C2 * (C3 / T / np.sinh(C3 / T)) ** 2 + C4 * (C5 / T / np.cosh(C5 / T)) ** 2)
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
        # Try to find ways to reduce computing time
        return H


class Stream:
    '''
        Describes material stream (IDEAL property method?)

        Attributes
        ----------
        - MW [kg/kgmol] Molar weight
        - RHOIG [kg/m3] Ideal gas density @ Actual conditions
        - STDRHOIG [kg/sm3] Ideal gas density @ Std. Conditions
        - CPIG [kJ/(mol*K)] Stream ideal gas heat capacity at constant pressure
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
        R = 8.31446261815324  # Gas Constant [J/(mole*K)]

        molweight = 0
        CPIG = 0
        massflow = 0
        for comp in self.compset:
            molweight += self.COMPMOLFR[comp.name] * comp.MW
            CPIG += self.COMPMOLFR[comp.name] * comp.CPIG(T)
            massflow += self.FLMOL * self.COMPMOLFR[comp.name] * comp.MW

        '''[kg/kgmol] Stream molar weight'''
        self.MW = molweight

        '''[kg/m3] Stream density @ Actual conditions'''
        self.RHOIG = self.MW / (R * self.T / (self.P * 1000))

        '''[kg/sm3] Stream density @ Standard conditions'''
        self.STDRHOIG = self.MW/ (R * self.T)

        '''[J/(mol*K)] Stream ideal gas heat capacity at constant pressure'''
        self.CPIG = CPIG

        '''[kg/hr] Stream mass flow'''
        self.FLMASS = massflow

        '''[m3/hr] Stream volume flow @ Actual conditions'''
        self.FLVOLIG = self.FLMOL * self.MW / self.RHOIG

        '''[sm3/hr] Stream volume flow @ Standard conditions'''
        self.STDVOLIG = self.FLMOL * self.MW/ self.STDRHOIG

        x_mass = dict()
        molconc = dict()
        indmolflows = dict()
        indmassflows = dict()
        for comp in self.compset:
            x_mass[comp.name] = self.COMPMOLFR[comp.name] * comp.MW / self.MW
            molconc[comp.name] = self.FLMOL * self.COMPMOLFR[comp.name] / self.FLVOLIG
            indmolflows[comp.name] = self.FLMOL * self.COMPMOLFR[comp.name]
            indmassflows[comp.name] = indmolflows[comp.name] * comp.MW

        '''[kgmol/hr] Molar flow of individual components'''
        self.FLINDMOL = indmolflows

        '''[mass. fract.] Stream composition in terms of mass fractions'''
        self.COMPMASSFR = x_mass

        '''[kgmol/m3] Stream composition in terms of molar concentrations @ Actual Conditions'''
        self.COMPMOLCIG = molconc

        '''[kg/hr] Mass flow of individual components'''
        self.FLINDMASS = indmassflows


class Reaction:
    '''
    Describes kinetic reactions

    Methods
    ----------
    .rate(T: float)
        Calculates Reaction Rate at specified temperature with Arrhenius Law
    '''
    def __init__(self, name: str, reagents: list[Species], stoic: list[float], dH: float, k0: float, E0: float):
        '''
        :param name: Reaction Name
        :param reagents: List of all reagents (same order as in equation)
        :param stoic: List of stoichiometric coefficients for listed components (same order as reagents). Negative values
        for reagents (left side of equation), positive for products (right side)
        :param dH: [kJ/mol] Heat of Reaction (Implemented temporarily before calculation as difference of enthalpies of formation
        would be introduced)
        :param k0: Reaction Rate Constant (Arrhenius Parameter)
        :param E0: [kJ/mol] Activation Energy
        '''
        self.name = name
        self.reagents = reagents
        compnameslist = list()
        for comp in reagents:
            compnameslist.append(comp.name)
        self.stoic = dict(zip(compnameslist, stoic))
        self.dH = dH
        self.k0 = k0
        self.E0 = E0

    def rate(self, T: float, conc: dict):
        '''
        Returns Reaction Rate at specified temperature

        :param T: [K] Temperature
        :param conc: [kmol/m3] Concentrations of components
        :return: Reaction Rate
        '''
        '''Multiplier that considers contribution of concentrations to Arrhenius Law'''
        mult = 1
        for comp in self.reagents:
            if self.stoic[comp.name] < 0:
                # Now calculates with concentrations in [kmol/m3]!
                mult = mult * ((conc[comp.name] / 1000) ** abs(self.stoic[comp.name]))
                '''Factor of 1000 added because initial units was mol/hr while
                concentrations calculated in kmol/hr'''
                #USE REACTION ORDER INSTEAD OF STOIC COEFFS !?
        # Needs to be revised!
        '''
        equation used: v = k * prod([I] ^ i)
            where k = k0 * exp(E0 / R / T)
        '''
        rate = mult * self.k0 * np.exp(-self.E0 / 8.3144 / T)
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
        :param numtubes: [] Number of reactor tubes
        :param rxnset: [list of Reaction] Set of reactions occurring in reactor
        '''
        self.length = length
        self.diameter = diameter
        self.numtubes = numtubes
        self.rxnset = rxnset

    def Simulation(self, inlet: Stream, dl: float, log: bool) -> tuple[Stream, pd.DataFrame]:
        '''
        Performs  integration along reactor length with Euler method for reactions listed

        :param inlet: [Stream] Reactor inlet stream
        :param step: [m] Integration resolution (step along reactor length)
        :param log: Select if tabular results for each timestep is required
        :return: [Stream] Reactor outlet stream and [pd.DataFrame] Calculations results on each iteration
        '''
        flow = inlet
        cell_volume = np.pi * (self.diameter ** 2) / 4 * dl  # [m3]
        l = 0  # [m]
        list_l = []
        list_T = []
        t = 0  # [s]
        temp_df = pd.DataFrame()
        while l < self.length:
            print('\tintegration l = {:.3f} m'.format(l))
            print('\t            t = {:.3f} s'.format(t))
            dt = cell_volume / flow.FLVOLIG * 3600  # [s]
            act_C = flow.COMPMOLCIG  # [kgmol/m3]
            act_T = flow.T  # [K]
            act_P = flow.P  # [MPa]
            act_Cp = flow.CPIG  # [kJ/(kg*K)]
            dT = 0
            ############################################
            '''matrices!!!'''
            ############################################
            generation = dict()
            for comp in flow.compset:
                generation[comp.name] = 0
            for rxn in self.rxnset:
                rate = rxn.rate(act_T, act_C)  # [kgmol/(m3*s)]
                dT += -rxn.dH * 1000 * rate  # [K]
                for comp in rxn.reagents:
                    generation[comp.name] += dt * rxn.stoic[comp.name] * rate * cell_volume  # [kgmol/s]
            new_indmolflows = dict()
            new_molfract = dict()
            new_molflow = flow.FLMOL + sum(generation.values())
            for comp in flow.compset:
                new_indmolflows[comp.name] = flow.FLINDMOL[comp.name] + generation[comp.name] * 3600 / 1000  # [kgmol/hr]
                new_molfract[comp.name] = new_indmolflows[comp.name] / new_molflow  # [mol. fract.]
            # new_T = act_T + dt * dT * 0.00845 * act_T / act_P / (act_Cp * 1000) # [K]
            new_T = act_T
            ### UNITS FOR HEAT BALANCE EQUATION?
            flow = Stream(flow.compset, new_molfract, new_molflow, act_P, new_T)
            temp_df = temp_df.append(flow.COMPMOLFR, ignore_index=True)
            list_l.append(l)
            list_T.append(new_T)
            l += dl
            t += dt
        temp_df['l'] = list_l
        temp_df['T'] = list_T
        temp_df = temp_df.set_index('l', drop=True)
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
    for param in results.columns[0:5]:
        # if 'x - ' in param:
        axes[0].plot(results.index, results[param], label= param)
    axes[0].legend()
    axes[1].set_ylabel('Temperature, [K]')
    axes[1].set_xlabel('Time, [s]')
    axes[1].plot(results.index, results['T'])
    plt.grid()
    plt.show()
    return 1



