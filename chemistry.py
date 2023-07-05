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
import math
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
        # print(self.name)
        for comp in self.reagents:
            if self.order[comp.ID] < 0:
                # print('comp ', comp.formula, 'conc=', conc[comp.ID])
                mult = mult * ((conc[comp.ID]) ** abs(self.order[comp.ID]))
                # print('mult=', mult)
                reactatnts_conc = np.append(reactatnts_conc, conc[comp.ID])
        # Needs to be revised (use matrix instead of loop)!
        max_rate = min(reactatnts_conc) / dt
        # print('mrate= ', max_rate)
        rate = mult * self.k0 * np.exp(-(self.E0 * 1000) / 8.3144 / T)  # [kgmol/(m3*s)]
        # print('calc rate= ', rate)
        if rate <= max_rate:
            pass
        else:
            # pass
            rate = max_rate
        # print('rate= ', rate)
        # print('consumed= ', rate * dt, end= '\n\n')


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
    def simulation(self, inlet: Stream, init_dl: float, log: bool) -> tuple[Stream, pd.DataFrame]:
        '''
        Performs  integration along reactor length with Euler method for reactions listed

        :param inlet: [Stream] Reactor inlet stream
        :param dl: [m] Integration resolution (step along reactor length)

        :param log: Select if tabular results for each timestep is required
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
            # print(rxn.reagents)
            # print(self.rxnset[7].reagents)
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
        # print(f'cell volume is {cell_volume}\ttube volume is {self.tubevolume}\t'
        #       f'rctr volume is {self.tubevolume * self.numtubes}\tcell duty is {cell_duty}')

        def step(cell_inlet: Stream, correction_factor: float, cell_dl: float, negative_balance_check_instep: bool) -> tuple[Stream, float, dict, int, bool]:
            cell_volume = np.pi * (self.diameter ** 2) / 4 * cell_dl# * self.numtubes # [m3]  ### ADD TUBES NO
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


            # '''Loop to update components concentrations after each rxn takes place'''
            # dQ = 0
            # rates_hist = dict()
            # for rxn in self.rxnset:
            #     '''Dictionary of stoichiometric coefficients'''
            #     stoic_dict = dict()
            #     for comp in cell_inlet.compset:
            #         if comp in rxn.reagents:
            #             stoic_dict[comp.ID] = rxn.stoic[comp.ID]
            #         else:
            #             stoic_dict[comp.ID] = 0
            #     '''Sort stoic dictionary and convert to vector'''
            #     stoic_vector = np.array(list(dict(sorted(stoic_dict.items())).values()))
            #     '''Calculate rxn rate'''
            #     rate = rxn.rate(act_T, act_C, dt)
            #     rates_hist[rxn.name] = rate
            #     '''Functional form of PFReactor mass balance differential equation for integration methods
            #     (It is declared each time inside the loop to allow its integration with m. methods without
            #     specifying any additional parameters)'''
            #     def concentrations_derivative_single(x, y, _stoic_vector=stoic_vector, _rate=rate, _T=act_T):
            #         x = 1
            #         return stoic_vector * rate
            #     '''Update components concentrations at cell outlet'''
            #     C_vect = m.integrate('rungekutta4th', concentrations_derivative_single, 1, C_vect, dt)
            #     act_C = dict(zip(comp_keys, C_vect))
            #     '''Add rxn heat to total heat at rctr cell'''
            #     dQ += rate * rxn.dH * -1000  # [kJ/(m3*s)]
            # # print(f'\tdQ at step {dQ}\t ext Q {cell_duty}', end= '')


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
                # NEGATIVE rxn RATES!
                return (_stoic_matrix * _rateconst_matrix).sum(axis= 0)

            '''Reactions rate constants matrix [No. rxns x 1] (for new concentrations?)'''
            rateconst_matrix = np.array(list(map(lambda x: x.rate(act_T, act_C, dt), self.rxnset)))  # [kgmol/(m3*s)]
            rates_hist = dict(zip(rxn_keys, rateconst_matrix))
            rateconst_matrix = np.reshape(rateconst_matrix, (len(rateconst_matrix), 1))

            '''Comps conc at cell outlet from PFReactor diff equation [1 x No. comps]'''
            # print(conc_frames)
            C_vect, cell_dl, integration_status = m.integrate('rungekuttamerson', concentrations_derivative, 1, C_vect, dt, negative_balance_check_instep)  # [kgmol/m3]
            # C_vect, cell_dl, integration_status = m.integrate('rungekutta4th', concentrations_derivative, 1, C_vect, dt, balance_check_instep)  # [kgmol/m3]
            # print('dl ->', dl)
            # C_vect, cell_dl, integration_status = m.integrate('gear', concentrations_derivative, step_frames, conc_frames, dt, balance_check_instep)  # [kgmol/m3]

            """
            act_C_aux = act_C
            rateconst_matrix1 = np.array(list(map(lambda x: x.rate(act_T, act_C_aux, dt), self.rxnset)))
            rateconst_matrix1 = np.reshape(rateconst_matrix1, (len(rateconst_matrix1), 1))
            k1 = (stoic_matrix * rateconst_matrix1).sum(axis= 0)
            i = 0
            for key in comp_keys:
                act_C_aux[key] = act_C_aux[key] + k1[i] * dt / 2
                i += 1
            rateconst_matrix2 = np.array(list(map(lambda x: x.rate(act_T, act_C_aux, dt), self.rxnset)))
            rateconst_matrix2 = np.reshape(rateconst_matrix2, (len(rateconst_matrix2), 1))
            k2 = (stoic_matrix * rateconst_matrix2).sum(axis= 0)
            i = 0
            for key in comp_keys:
                act_C_aux[key] = act_C_aux[key] + k2[i] * dt / 2
                i += 1
            rateconst_matrix3 = np.array(list(map(lambda x: x.rate(act_T, act_C_aux, dt), self.rxnset)))
            rateconst_matrix3 = np.reshape(rateconst_matrix3, (len(rateconst_matrix3), 1))
            k3 = (stoic_matrix * rateconst_matrix3).sum(axis= 0)
            i = 0
            for key in comp_keys:
                act_C_aux[key] = act_C_aux[key] + k2[i] * dt
                i += 1
            rateconst_matrix4 = np.array(list(map(lambda x: x.rate(act_T, act_C_aux, dt), self.rxnset)))
            rateconst_matrix4 = np.reshape(rateconst_matrix4, (len(rateconst_matrix4), 1))
            k4 = (stoic_matrix * rateconst_matrix4).sum(axis= 0)
            C_vect = C_vect + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            """

            # Different results using m.integrate() and explicitly writing 4th order method!

            '''Update comps concentration dictionary'''
            act_C = dict(zip(comp_keys, C_vect))
            # print(max(act_C.values()))
            # print('\n', act_C, '\n')
            # print()
            # print(act_C)

            '''Sum of reaction heat for all rxns in rctr [1 x 1]'''
            dQ = np.sum(rateconst_matrix * rxndH_matrix) * -1000  # [kJ/(m3*s)]

            '''Functional form of PFReactor heat balance differential equation for integration methods'''
            def temperature_derivative(x, y, _dQ=dQ, _P=act_P, _Cp=act_Cp):
                x = 1
                return (_dQ + cell_duty) * 0.008314463 * y / _P / _Cp
            '''Update cell temperature'''
            # new_T = m.integrate('rungekutta4th', temperature_derivative, 1, act_T, dt, balance_check_instep)
            new_T, cell_dl_temperature, integration_status_temperature = m.integrate('rungekuttamerson', temperature_derivative, 1, act_T, dt, negative_balance_check_instep)

            # new_T = act_T
            '''Comps mole fractions list to replace negative values'''
            # print('\n', act_C)
            act_molfract = list(map(lambda x: act_C[x] / sum(act_C.values()), comp_keys))
            # act_molfract = list([num if num > 0 else 0 for num in act_molfract])
            '''Comps mole fractions at cell outlet'''
            new_compmolfr = dict(zip(comp_keys, act_molfract))  # [mol. fract.]
            '''Comps mole flow at cell outlet (volume calculated from PR EOS or IG EOS at cell inlet)'''
            new_molflow = sum(list(map(lambda x: volflow * act_C[x], comp_keys)))  # [kgmol/hr]
            '''Update flow to get estimated flow at cell outlet'''
            outlet = Stream(flow.compset, new_compmolfr, new_molflow, act_P, new_T, inlet.eos_option)
            # print('\t', sum(cell_outlet.COMPMASSFR.values()))
            if any(x < 0 for x in outlet.COMPMASSFR.values()):
                # print('\nERROR - Negative mass fractions obtained!')
                outlet = cell_inlet
                integration_status = 1
                negative_balance_check_instep = True
                # sys.exit()

            return outlet, dt, rates_hist, integration_status, negative_balance_check_instep



        '''Integration through reactor length'''
        dl = init_dl
        while l <= self.length:
            # Printout progress bar
            # sys.stdout.write('\tintegrating at {:.4e}/{:.2f} m with step of {:.2e} m\n\r'.format(l, self.length, dl))
            sys.stdout.write('\tintegrating at {:.4f}/{:.4f} m with step of {:.2e} m\r'.format(l, self.length, dl))
            negative_balance_check = False
            counter = 0
            step_status = -1
            while not step_status == 0:
                # sys.stdout.write('balance status: {}\r'.format(balance_check))
                cell_outlet, dt, rates_hist, step_status, negative_balance_check = step(flow, 1, dl, negative_balance_check)
                # print(max(cell_outlet.COMPMOLC.values()), dl)
                # print(step_status, dl)


                if step_status == 1:
                    dl *= 0.8
                    counter += 1
                    # print('whoops')
                elif step_status == 2:
                    dl *= 1.2
                # print(counter)
                # if l >= 3.e-7:
                #     dl *= 1.01
                # if 6e-8 <= l <= 6.1e-8:
                #     dl *= 0.1
                # if 6.3e-8 <= l <= 6.31e-8:
                #     dl *= 10


                ###!!!########################################################
                mass_balance_check = abs(flow.FLMASS - cell_outlet.FLMASS) <= 1e-5
                # mass_balance_check = True
                corrector_l = 0.5
                corrector_r = 2
                cell_outlet_m = cell_outlet
                while not mass_balance_check:
                    corrector_m = (corrector_l + corrector_r) / 2
                    cell_outlet_l = Stream(cell_outlet.compset, cell_outlet.COMPMOLFR, cell_outlet.FLMOL * corrector_l,
                                           cell_outlet.P, cell_outlet.T, 'IG')
                    cell_outlet_m = Stream(cell_outlet.compset, cell_outlet.COMPMOLFR, cell_outlet.FLMOL * corrector_m,
                                           cell_outlet.P, cell_outlet.T, 'IG')
                    cell_outlet_r = Stream(cell_outlet.compset, cell_outlet.COMPMOLFR, cell_outlet.FLMOL * corrector_r,
                                           cell_outlet.P, cell_outlet.T, 'IG')
                    check_l = cell_outlet_l.FLMASS - flow.FLMASS
                    check_m = cell_outlet_m.FLMASS - flow.FLMASS
                    check_r = cell_outlet_r.FLMASS - flow.FLMASS
                    if abs(check_m) <= 1e-5:
                        mass_balance_check = True
                    else:
                        if check_m * check_l < 0:
                            corrector_r = corrector_m
                        elif check_m * check_r < 0:
                            corrector_l = corrector_m
                        else:
                            print('Here we have a problem!')
                cell_outlet = cell_outlet_m
                if not mass_balance_check:
                    print(mass_balance_check)
                ###!!!##########################################################

            """
            '''Setup variables to control material balance convergence loop'''
            correction_l = 1 * 0.5
            correction_r = 1 * 1.5
            mbal = False
            counter = 0
            # Bypass mbal check
            mbal = True
            cell_outlet_m, dt_m, rate_hist_m = step(flow, 1)
            while not mbal:
                counter += 1
                correction_m = (correction_l + correction_r) / 2
                cell_outlet_l, dt_l, rate_hist_l = step(flow, correction_l)
                cell_outlet_m, dt_m, rate_hist_m = step(flow, correction_m)
                cell_outlet_r, dt_r, rate_hist_r = step(flow, correction_r)
                check_l = cell_outlet_l.FLMASS - flow.FLMASS
                check_m = cell_outlet_m.FLMASS - flow.FLMASS
                check_r = cell_outlet_r.FLMASS - flow.FLMASS
                if abs(check_m) <= 1e-4:
                    mbal = True
                else:
                    if check_m * check_l < 0:
                        correction_r = correction_m
                    elif check_m * check_r < 0:
                        correction_l = correction_m
                    else:
                        print('Here we have a problem!')
            # print(f'\tMBAL converged in {counter} iterations')
            cell_outlet, dt, rates_hist = cell_outlet_m, dt_m, rate_hist_m
            """

            '''Update flow to cell outlet conditions'''
            flow = cell_outlet
            '''Step forward through reactor'''
            l += dl
            t += dt
            '''Dict to store output variables inside the loop before appending to results dataframe'''
            output_line = dict()
            '''Fill output_line'''
            # comp_keys_output= list(map(lambda x: x, flow.compset))
            comp_conc_output = flow.COMPMOLFR.copy().values()
            output_line.update(dict(zip(comp_keys_output, comp_conc_output)))
            output_line['FLMOL [kgmol/hr]'] = flow.FLMOL
            output_line['l [m]'] = l
            output_line['t [s]'] = t
            output_line['T [K]'] = flow.T
            # output_line.update(rates_hist)
            '''Add output line as new index to list of frames'''
            frames.append(pd.DataFrame([output_line]))
            conc_frame = np.array(list(flow.COMPMOLC.copy().values()))
            step_frame = l
            conc_frames = np.append(conc_frames, [conc_frame], axis= 0)

            step_frames = np.append(step_frames, step_frame)
        '''Combine all frames to DataFrame'''
        output_df = pd.concat(frames)
        '''Set lengths column as df indexes for result df'''
        output_df = output_df.set_index('l [m]', drop=True)

        print(f'Reactor {"INSERTYOURREACTORIDHERE"} integration completed:\n'
              f'\trctr length {self.length : .2f} m\t\t\tresidense time {t * 1000 : .2f} ms\n'
              f'\tlast length step {dl : .3e} m\t\tlast time step {dt * 1000 : .3e} ms')

        return flow, output_df
