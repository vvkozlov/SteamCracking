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
        [5] Aspen HYSYS User Manual

"""
import math
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import usermath
from coreobjects import Species, Stream
import usermath as m
from progress.bar import IncrementalBar


class Reaction:
    """
    Describes kinetic reactions

    Methods
    ----------
    .rate(T: float)
        Calculates Reaction Rate at specified temperature with Arrhenius Law
    """
    def __init__(self, ID:int, name: str, reagents: list[Species], stoic: list[float], order: list[float],
                 dH: float, k0: float, E0: float, sequence: int):
        """
        :param ID: Reaction ID
        :param name: Reaction Name
        :param reagents: List of all reagents (same order as in equation)
        :param stoic: List of stoichiometric coefficients for listed components (written in same order as reagents).
        Negative values for reagents (left side of equation), positive for products (right side)
        :param order: List of reaction orders by each reactant for listed components (written in same order as reagents).
        Negative values for reagents (left side of equation), positive for products (right side)
        :param dH: [kJ/mol] Heat of Reaction (if equals zero value will be obtained from reagents enthalpy difference)
        :param k0: Reaction Rate Constant (Arrhenius Parameter)
        :param E0: [kJ/mol] Activation Energy
        :param sequence: Reaction type (1 - initiation, 2 - propagation, 3 - termination) (TEST MODE)
        """
        self.ID = ID
        self.name = name
        self.reagents = reagents
        self.stoic = dict(zip(list(map(lambda x: x.ID, self.reagents)), stoic))
        self.order = dict(zip(list(map(lambda x: x.ID, self.reagents)), order))
        self.k0 = k0
        self.E0 = E0
        self.sequence = sequence

        DHFORM_vect = np.array(list(map(lambda x: x.DHFORM, reagents)))
        if dH == 0:
            self.dH = np.sum(np.array(stoic) * DHFORM_vect) / 10**6  # [kJ/mol]
        else:
            self.dH = dH

    def rate(self, T: float, conc: dict, dt: float):
        '''
        Returns Reaction Rate for forward rxns at specified temperature

        :param T: [K] Temperature
        :param conc: [kmol/m3] Concentrations of components
        :return: [kgmol/(m3*s)] Reaction Rate
        '''
        '''
        equation used: v = k * prod([I] ^ i)
            where k = k0 * exp(E0 / R / T)
        '''
        mult = 1  # [kgmol/m3]^n Multiplier that considers contribution of concentrations
        reactatnts_conc = np.array([])
        for comp in self.reagents:
            if self.order[comp.ID] < 0:
                mult = mult * ((conc[comp.ID]) ** abs(self.order[comp.ID]))
                reactatnts_conc = np.append(reactatnts_conc, conc[comp.ID])
        # Needs to be revised (use matrix instead of loop)!
        max_rate = min(reactatnts_conc) / dt
        rate = mult * self.k0 * np.exp(-(self.E0 * 1000) / 8.3144 / T)  # [kgmol/(m3*s)]
        if rate <= max_rate:
            pass
        else:
            rate = max_rate



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
    def __init__(self, length: float, diameter: float, numtubes: float, rxnset: list[Reaction],
                 duty: float = 0):
        '''
        :param length: [m] Reactor tube length
        :param diameter: [m] Reactor tube diameter
        :param numtubes: [No] Number of reactor tubes
        :param rxnset: [list of Reaction] Set of reactions occurring in reactor
        :param duty: [kW] External rctr heat duty (furnace heat flow)
        '''
        self.length = length
        self.diameter = diameter
        self.numtubes = numtubes
        self.tubevolume = np.pi * (self.diameter ** 2) / 4 * self.length
        self.rxnset = rxnset
        self.duty = duty

    # without matrices 1863.49 ms
    def simulation(self, inlet: Stream, init_dl: float, log: bool, output_reduction_factor: int) -> tuple[Stream, pd.DataFrame]:
        '''
        Performs  integration along reactor length with Euler method for reactions listed

        :param inlet: [Stream] Reactor inlet stream
        :param dl: [m] Integration resolution (step along reactor length)
        :param log: Select if tabular results for each timestep is required
        :param output_reduction_factor: Reduces number of output lines in logfine by factor of output_reduction_factor
        :return: [Stream] Reactor outlet stream and [pd.DataFrame] Calculations results on each iteration
        '''
        '''Determine conditions at rctr inlet'''
        dl = init_dl
        flow = inlet
        l = 0  # [m]
        t = 0  # [s]
        '''Store convergence limit for mass balance iterations'''
        mbal_tol = 1e-3
        '''External heat duty per rctr cell considering equal heat distribution through rctr volume'''
        cell_duty = self.duty / (self.tubevolume * self.numtubes)  # [kJ/(m3*s)]
        '''Keys for components and reactions lists - to make sure that all matrices are uniform'''
        comp_keys = sorted(flow.COMPMOLFR)  # Some more elegant way to create matching list should be found
        aux_comp_list = list(map(lambda x: x, flow.compset))
        aux_comp_list.sort(key= lambda x: x.ID, reverse= False)
        comp_keys_output = list(map(lambda x: x.formula, aux_comp_list))
        rxn_keys = list(map(lambda x: x.ID, self.rxnset))
        rxn_keys_output = list(map(lambda x: x.name, self.rxnset))
        '''Reaction stoich coefficients matrix [No. rxns x No. comps]'''
        stoic_df = pd.DataFrame(index=rxn_keys, columns=comp_keys)
        '''Reaction order matrix [No. rxns x No. comps]'''
        order_df = pd.DataFrame(index=rxn_keys, columns=comp_keys)
        '''Reaction enthalpy difference dH matrix [No. rxns x 1]'''
        rxndH_df = pd.DataFrame(index=rxn_keys, columns=['dH'])
        '''Assemble stoich coeffs and rxn enthalpies df's'''
        for rxn in self.rxnset:
            rxndH_df['dH'][rxn.ID] = rxn.dH
            for comp in inlet.compset:
                if comp in rxn.reagents:
                    stoic_df[comp.ID][rxn.ID] = rxn.stoic[comp.ID]
                else:
                    stoic_df[comp.ID][rxn.ID] = 0
        '''Convert stoich coeffs and rxn enthalpies df's to matrices'''
        stoic_matrix = np.array(stoic_df)
        rxndH_matrix = np.array(rxndH_df)
        '''Create storages for frames of output DataFrame '''
        frames = []
        conc_frames = np.array([np.array(list(dict(sorted(inlet.COMPMOLC.items())).values()))])
        step_frames = np.array([0])
        def step(cell_inlet: Stream, correction_factor: float, cell_dl: float, negative_balance_check_instep: bool) -> tuple[Stream, float, dict, int, bool]:
            cell_volume = np.pi * (self.diameter ** 2) / 4 * cell_dl #* self.numtubes # [m3]  ### ADD TUBES NO
            '''Calculate volume flowrate trough cell and determine initial concentrations'''
            volflow = cell_inlet.FLVOL * correction_factor
            act_C = cell_inlet.COMPMOLC  # [kgmol/m3]
            '''Residence time for finite rctr cell'''
            dt = cell_volume / volflow * 3600  # [s]
            '''Determine conditions at cell inlet'''
            act_T = cell_inlet.T  # [K]
            act_P = cell_inlet.P  # [MPa]
            act_Cp = cell_inlet.CP  # [J/(mol*K)] - only IG option availible
            '''Comps concentrations vector [1 x No. comps]'''
            C_vect = np.array(list(dict(sorted(act_C.items())).values()))  # [kgmol/m3]

            '''Functional form of PFReactor mass balance differential equation for integration methods'''
            def concentrations_derivative(x, y, _dt=None, _stoic_matrix=None, _rxnset=None, _T=None):
                if _rxnset is None:
                    _rxnset = self.rxnset
                if _stoic_matrix is None:
                    _stoic_matrix = stoic_matrix
                if _T is None:
                    _T = act_T
                if _dt is None:
                    _dt = dt
                x = 1
                '''Reactions rate constants matrix [No. rxns x 1]'''
                _rateconst_matrix = np.array(list(map(lambda x: x.rate(_T, dict(zip(comp_keys, y)), _dt), _rxnset)))  # [kgmol/(m3*s)]
                _rateconst_matrix = np.reshape(_rateconst_matrix, (len(_rateconst_matrix), 1))
                return (_stoic_matrix * _rateconst_matrix).sum(axis= 0)

            '''Reactions rate constants matrix [No. rxns x 1] (for new concentrations?)'''
            rateconst_matrix = np.array(list(map(lambda x: x.rate(act_T, act_C, dt), self.rxnset)))  # [kgmol/(m3*s)]
            rates_hist = dict(zip(rxn_keys, rateconst_matrix))
            rateconst_matrix = np.reshape(rateconst_matrix, (len(rateconst_matrix), 1))

            '''Comps conc at cell outlet from PFReactor diff equation [1 x No. comps]'''
            C_vect, cell_dl, integration_status, LTE_instep = m.integrate('rungekuttamerson_adaptive', concentrations_derivative, 1, C_vect, dt, negative_balance_check_instep)  # [kgmol/m3]

            # negative_balance_check_instep = True
            # C_vect, cell_dl, integration_status, LTE_instep = m.integrate('rungekuttamerson_fixed', concentrations_derivative, 1, C_vect, dt, negative_balance_check_instep)  # [kgmol/m3]
            # ADDITIONAL VARIABLE 'LTE' ADDED TO SCREEN LTE
            # print('dl ->', dl)
            # C_vect, cell_dl, integration_status = m.integrate('gear', concentrations_derivative, step_frames, conc_frames, dt, negative_balance_check_instep)  # [kgmol/m3]
            # Different results using m.integrate() and explicitly writing 4th order method!

            '''Update comps concentration dictionary'''
            act_C = dict(zip(comp_keys, C_vect))

            '''Sum of reaction heat for all rxns in rctr [1 x 1]'''
            dQ = np.sum(rateconst_matrix * rxndH_matrix) * -1000  # [kJ/(m3*s)]

            '''Functional form of PFReactor heat balance differential equation for integration methods'''
            def temperature_derivative(x, y, _dQ=dQ, _P=act_P, _Cp=act_Cp):
                x = 1
                return (_dQ + cell_duty) * 0.008314463 * y / _P / _Cp
            '''Update cell temperature'''
            # new_T = m.integrate('rungekutta4th', temperature_derivative, 1, act_T, dt, balance_check_instep)
            new_T, cell_dl_temperature, integration_status_temperature, LTE_temperature = m.integrate('rungekuttamerson_adaptive', temperature_derivative, 1, act_T, dt, negative_balance_check_instep)
            # new_T, cell_dl_temperature, integration_status_temperature, LTE_temperature = m.integrate('rungekuttamerson_fixed', temperature_derivative, 1, act_T, dt, negative_balance_check_instep)
            # new_T = act_T
            '''Comps mole fractions list to replace negative values'''
            act_molfract = list(map(lambda x: act_C[x] / sum(act_C.values()), comp_keys))
            '''Comps mole fractions at cell outlet'''
            new_compmolfr = dict(zip(comp_keys, act_molfract))  # [mol. fract.]
            '''Comps mole flow at cell outlet (volume calculated from PR EOS or IG EOS at cell inlet)'''
            new_molflow = sum(list(map(lambda x: volflow * act_C[x], comp_keys)))  # [kgmol/hr]
            '''Update flow to get estimated flow at cell outlet'''
            outlet = Stream(flow.compset, new_compmolfr, new_molflow, act_P, new_T, inlet.eos_option)
            # print('\t', sum(cell_outlet.COMPMASSFR.values()))
            if any(x < 0 for x in outlet.COMPMASSFR.values()):
                # print('ERROR - Negative mass fractions obtained!\t\t\t\t', end = '\r')
                outlet = cell_inlet
                integration_status = 1
                negative_balance_check_instep = True
                # sys.exit()

            return outlet, dt, rates_hist, integration_status, negative_balance_check_instep, LTE_instep



        '''Integration through reactor length'''
        dl = init_dl
        steps_counter = 0
        while l <= self.length:
            '''Printout progress bar'''
            sys.stdout.write('\tintegrating at {:.4f}/{:.2f} m with step of {:.2e} m\r'.format(l, self.length, dl))
            '''Check for reactions balance (if calculated comps consumption rate exceeds actual flowrate in cell)'''
            negative_balance_check = False
            ode_conv_counter = 0
            step_status = -1
            while not step_status == 0:
                # sys.stdout.write('balance status: {}\r'.format(balance_check))
                cell_outlet, dt, rates_hist, step_status, negative_balance_check, LTE = step(flow, 1, dl, negative_balance_check)
                if step_status == 1:
                    dl *= 0.8
                elif step_status == 2:
                    dl *= 1.2
                ode_conv_counter += 1
                '''Tune volume flow to converge mass balance'''
                cell_outlet, mbal_conv_counter = usermath.mbal_check('bisection', flow, cell_outlet, 1e-5)
            steps_counter += 1

            '''Update flow to cell outlet conditions'''
            flow = cell_outlet
            '''Step forward through reactor'''
            l += dl
            t += dt
            '''Dict to store output variables inside the loop before appending to results dataframe'''
            output_line = dict()
            '''Fill output_line'''
            output_line['Length step'] = steps_counter
            output_line['ODE convergence steps'] = ode_conv_counter
            output_line['MBAL convergence steps'] = mbal_conv_counter
            comp_conc_output = flow.COMPMOLFR.copy().values()
            output_line.update(dict(zip(comp_keys_output, comp_conc_output)))
            output_line['FLMOL [kgmol/hr]'] = flow.FLMOL
            output_line['l [m]'] = l
            output_line['t [s]'] = t
            output_line['T [K]'] = flow.T
            output_line['LTE'] = LTE
            rates_hist_output = rates_hist.copy().values()
            output_line.update(dict(zip(rxn_keys_output, rates_hist_output)))
            '''Add output line as new index to list of frames'''
            frames.append(pd.DataFrame([output_line]))
            conc_frame = np.array(list(flow.COMPMOLC.copy().values()))
            step_frame = l
            conc_frames = np.append(conc_frames, [conc_frame], axis= 0)

            step_frames = np.append(step_frames, step_frame)
        '''Combine all frames to DataFrame'''
        output_df = pd.concat(frames)
        '''Reduce output df by factor of output_reduction_factor'''
        output_df = output_df[output_df['Length step'] % output_reduction_factor == 0]
        '''Set lengths column as df indexes for result df'''
        output_df = output_df.set_index('l [m]', drop=True)

        print(f'Reactor {"INSERTYOURREACTORIDHERE"} integration completed:\n'
              f'\trctr length {self.length : .2f} m\t\t\tresidense time {t * 1000 : .2f} ms\n'
              f'\tlast length step {dl : .3e} m\t\tlast time step {dt * 1000 : .3e} ms')
        return flow, output_df
