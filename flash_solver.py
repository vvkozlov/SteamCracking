import os
import eos_interfaces_v1 as intrf
import eos_engine_v4 as calc
import pandas as pd
import numpy as np
import time
import sys
import database_v2 as intdb
from rctrs_engine_v3 import Stream

'''
REFERENCES
- A. Jebarjadi - Multi-Phase Multi-Component Equilibrium Flash Calculations for CompFlow Bio using Modified Volume-Translated Peng-Robinson EOS (2017)
- GPSA Engineering Data Book 13th Edition
- D. Peng, D. Robinson - A New Two-Constant Equation of State
- M. Abramowitz - Handbook of Mathematical Functions with Formulas, Graphs and Matematical Tables
- HYSYS Components Properties Database
- Densities - Costald (Source?)
'''

'''Intermediate file to match formates of streams in rctr models and eos_engine'''


### Gas Constant
# R_SI = 8.31446261815324 # [J/(mole*K)] - SI
# R_field = 10.731577089016 # [psi*ft3/(lbmol*R)] - Field
convcrit = 10**-4 # Convergence criteria for K-values
steps_limit = 25


### Import of Components Properties and Binary Interraction Coefficients Databases
cwd = os.getcwd()
comppropDB = pd.read_excel(intrf.get_comppropDB_names(cwd)[0], index_col= 'Name')
binarycoefDB = pd.read_excel(intrf.get_binarycoefDB_names(cwd)[0], index_col= 'index_col')


### Input of main parameters
streamcomp_names_list = intrf.get_streamcomp_names(cwd)
print('Select stream composition file:')
streamcomp_name = 'Stream1 (GPSA exmpl) - StreamComposition.xlsx' # intrf.select_file(streamcomp_names_list) #
input_streamcomp = pd.read_excel(streamcomp_name, index_col= 'Name')

CL2 = intdb.CL2
PROPYLENE = intdb.PROPYLENE
C3H5CL = intdb.C3H5CL
HCL = intdb.HCL
C3H6CL2 = intdb.C3H6CL2
compset = [CL2, PROPYLENE, C3H5CL, HCL, C3H6CL2]
comp_x0 = dict({CL2.name : 0.2, PROPYLENE.name : 0.8, C3H5CL.name : 0, HCL.name : 0, C3H6CL2.name: 0})
molflow = 0.627017008965004  # Feed stream molar flow [kgmol/hr]
P = 0.2  # Reaction Pressure [MPa]
T0 = 200 + 273.15  # Initial Temperature [K]
inlet = Stream(compset, comp_x0, molflow, P, T0)

indexes = list()
for comp in inlet.compset:
	indexes.append(comp.name)
input_streamcomp = pd.DataFrame(index= indexes, columns= ['Content [mol. fract.]'])
for comp in inlet.compset:
	input_streamcomp.loc[comp.name]['Content [mol. fract.]'] = inlet.COMPMOLFR[comp.name]


print('Input operating pressure P [bara]:')
P = calc.UnitsConverter.Pressure.bar_to_kPa(float(input()))
P_field = calc.UnitsConverter.Pressure.kPa_to_psi(P)
print('Input operating temperature T [C]:')
T = float(input())
T_field = calc.UnitsConverter.Temperature.C_to_R(T)

print('\nINPUT SUMMARY:')
print('\t- Pressure - {:.1f} kPa\n\t- Temperature - {:.1f} C'.format(P, T))
print('\t- Stream composition [mol. fract.]:')
for component in input_streamcomp.index:
	print('\t\t{:>10} - {:.4f}'.format(component, input_streamcomp.loc[component]['Content [mol. fract.]']))
print('-'*30, end= '\n\n')

# print('Confirm input? ("y" - yes, "n" - no)')
# confirm = str(input())
# if confirm == 'n':
# 	sys.exit()

start_time = time.perf_counter()
pd.options.display.float_format = '{:>8.4f}'.format

### Performing Flash Calculations
# print('-->Input comp.:', input_streamcomp)
'''equicomp_df, phase_fractions, zfactors, K_values = calc.flash_calc_PR_EOS(comppropDB,
																binarycoefDB,
																input_streamcomp,
																P_field,
																T_field,
																convcrit,
																steps_limit)'''
z = calc.flash_calc_PR_EOS_z_only(inlet,
								  comppropDB,
								  binarycoefDB,
								  convcrit,
								  steps_limit)
print(z)
print('-->Output comp.:', equicomp_df)
print('-->Output phase fract.:', phase_fractions)
print('-->Output zfactors:', zfactors)
'''
print('\n\nPhase Compositions:\n', equicomp_df)
### Check phase compositions sums
composition_check = pd.DataFrame({(equicomp_df['vapor'].sum(),
								   equicomp_df['liquid'].sum(),
								   equicomp_df['aqueous'].sum())},
								 columns= (['vapor', 'liquid', 'aqueous']),
								 index= ['TOTAL'])

print(composition_check)
print('\nVapor fraction (mol.): {:.5f}\nLiquid Fraction (mol.): {:.5f}\nAqueous Fraction (mol.): {:.5f}'.format(phase_fractions['vapor'],
																											phase_fractions['liquid'],
																											phase_fractions['aqueous']))

'''
'''
### PHYSICAL PROPERTIES CALCULATIONS
print('\n\nPhase properties:')

vapor_density = calc.get_vapor_phase_density(comppropDB, equicomp_df, P_field, T_field, zfactors.loc['vapor']['zj'])
liquid_density = calc.get_liquid_phase_density(comppropDB, equicomp_df, T_field, 'liquid')

if input_streamcomp.loc['H2O']['Content [mol. fract.]'] == 0:
	phases_num = 2
	aqueous_density = np.nan
else:
	phases_num = 3
	aqueous_density = calc.get_liquid_phase_density(comppropDB, equicomp_df, T_field, 'aqueous')
'''
'''
densities = calc.get_phase_densities_actcond(comppropDB,
											 equicomp_df,
											 phase_fractions,
											 P_field,
											 T_field,
											 zfactors)
print('Vapor phase Density [kg/m3]:\t {:.2f}'.format(densities['vapor']))
print('Liquid phase Density [kg/m3]:\t {:.2f}'.format(densities['liquid']))
print('Aqueous phase Density [kg/m3]:\t {:.2f}'.format(densities['aqueous']))

phaseMW = pd.Series([calc.get_phase_molar_weigh(comppropDB, equicomp_df, 'vapor'),
					 calc.get_phase_molar_weigh(comppropDB, equicomp_df, 'liquid'),
					 calc.get_phase_molar_weigh(comppropDB, equicomp_df, 'aqueous')], index=['V', 'L', 'Q'])
# print('mixro={}'.format(calc.get_mix_density(densities, phase_fractions, phaseMW)))

### DOES NOT WORK STABLE WITH WATER AD DOES NOT WORK WITHOUT WATER !!!

# print('MW vapor:\t{}\nMW liquid:\t{}\nMW aqueous:\t'.format(calc.get_phase_molar_weigh(comppropDB, equicomp_df, 'vapor'),
# 															calc.get_phase_molar_weigh(comppropDB, equicomp_df, 'liquid')))
'''

print('\nExecution time: {:.1f} ms\n'.format((time.perf_counter() - start_time) * 1000))

'''
try:
	writer = pd.ExcelWriter('solver_output.xlsx')
	equicomp_df.to_excel(writer)
	writer.save()
except:
	print('Failed to write down results')


print('Open output file? ("y" - yes, "n" - no)')
confirm = str(input())
if confirm == 'y':
	try:
		os.system('start excel.exe "%s//solver_output.xlsx"' % (sys.path[0], ))
	except:
		print('Failed to open file')
'''

check_df = ([0.919783235, 0.047922284, 0.014828746, 0.003731263, 0.002115643, 0.000428368, 0.000361431, 0.000182001, 4.52983E-05, 0, 0, 0, 0.010600829, 0])
print('\nTEST PASSED: {}'.format((abs(equicomp_df['vapor'] - check_df) < 10e-6).all()))