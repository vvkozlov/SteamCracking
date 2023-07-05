"""
Header      : usermath.py
Created     : 08.01.2023
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Contains some custom functions

References	:
			[1] - Add references for integration methods

"""
import sys
import pandas as pd
import numpy as np

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


"""Calculation of predictor and corrector coefficients for Gear implicit method"""
ns = 6  # Order of method
r = ns + 1
a_corrector = pd.DataFrame()
for s in range(2, r+1):
	ai = []
	for i in range(0, s):
		ai.append(0)
		for k in range(0, s):
			if k != i:
				prod = 1
				for j in range(0, s):
					if j != i and j != k:
						prod *= (s - j - 1) / (i - j)
					else:
						pass
				ai[i] += 1 / (i - k) * prod
			else:
				pass
	frame = pd.DataFrame({s-1: ai}).transpose()
	a_corrector = a_corrector.append(frame)
a_corrector = np.array(a_corrector.loc[ns])  # Vector of corrector coefficients of ns order
a_corrector = a_corrector[:, np.newaxis]


a_predictor = pd.DataFrame()
for s in range(2 , r+1):
	ai = []
	for i in range(0, s):
		ai.append(0)
		for k in range(0, s):
			if k != i:
				prod = 1
				for j in range(0, s):
					if j != i and j != k:
						prod *= (s - j) / (i - j)
					else:
						pass
				ai[i] += 1 / (i - k) * prod
			else:
				pass
	frame = pd.DataFrame({s-1: ai}).transpose()
	a_predictor = a_predictor.append(frame)
a_predictor = np.array(a_predictor.loc[ns])  # Vector of predictor coefficients of ns order
a_predictor = a_predictor[:, np.newaxis]
print('I am called!')  # Check how many times these vectors are calculated


def integrate(method: str, f, x0, y0, h: float, balance_check_inintegrator: bool) -> tuple:
	"""
		Calculate one step forward numerical integration using different methods.
		Available methods:
			- Euler's ('euler')
			- Runge-Kutta 4th order ('rungekutta4th')
			- Runge-Kutta-Felberg 5th order ('rungekuttafelberg5th')
		-----
		:param method: Key to select integration method
		:param f: Function-like expression to integrate through
		:param x0: X coordinate of initial point
		:param y0: Y coordinate of initial point
		:param h: Integration step
		"""
	if method == 'euler':
		return y0 + h * f(x0, y0)
	elif method == 'rungekutta4th':
		k1 = f(x0, y0)
		k2 = f(x0 + h / 2, y0 + k1 * h / 2)
		k3 = f(x0 + h / 2, y0 + k2 * h / 2)
		k4 = f(x0 + h, y0 + k3 * h)
		return y0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
	elif method == 'rungekuttafelberg5th':
		k1 = h * f(x0, y0)
		k2 = h * f(x0 + 1 / 4 * h, y0 + 1 / 4 * k1)
		k3 = h * f(x0 + 3 / 8 * h, y0 + 3 / 32 * k1 + 9 / 32 * k2)
		k4 = h * f(x0 + 12 / 13 * h, y0 + 1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3)
		k5 = h * f(x0 + h, y0 + 439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4)
		k6 = h * f(x0 + 1 / 2 * h, y0 - 8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5)
		return y0 + 16 / 135 * k1 + 6656 / 12825 * k3 + 28561 / 56430 * k4 - 9 / 50 * k5 + 2 / 55 * k6
	elif method == 'rungekuttamerson':
		'''termination_status stores the reason of integration termination:
			0 - good accuracy
			1 - error is too large, reduce dh
			2 - error is too small, increase dh'''
		termination_status = 0
		# !!! Feedback to 'dt' in simulation required for correct work !!!
		tolerance = 1e-11
		check = False
		h0 = h
		while not check:
			k0 = h0 * f(x0, y0)
			k1 = h0 * f(x0 + 1 / 3 * h0, y0 + 1 / 3 * k0)
			k2 = h0 * f(x0 + 1 / 3 * h0, y0 + 1 / 6 * k0 + 1 / 6 * k1)
			k3 = h0 * f(x0 + 1 / 2 * h0, y0 + 1 / 8 * k0 + 3 / 8 * k2)
			k4 = h0 * f(x0 + h0, y0 + 1 / 2 * k2 - 3 / 2 * k2 + 2 * k3)
			y = y0 + (k0 + 4 * k3 + k4) / 6
			R = (2 * k0 - 9 * k2 + 8 * k3 - k4) / 30
			try:
				R = max(abs(R))
			except:
				R = abs(R)
			# print(R)
			check1 = R <= tolerance
			check2 = R >= (tolerance / 30)
			'''Check if step status is changed by mass balances'''
			check3 = balance_check_inintegrator
			# print(max(abs(R)))
			# print('check1 - {}\tcheck2 - {}'.format(check1, check2))
			if check1 and check2:
				break
			elif not check1:
				termination_status = 1
				return y, h0 * 0.75, termination_status
			# use h0 * 0.5
			elif not check2:
				if not check3:
					termination_status = 2
					return y, h0 * 1.5, termination_status
				else:
					break
		return y, h0, termination_status
	elif method == 'gear':
		tolerance1 = 1e-7
		if len(x0) < 6:
			# k1 = h * f(x0, y0[-1])
			# k2 = h * f(x0[-1] + 1 / 4 * h, y0[-1] + 1 / 4 * k1)
			# k3 = h * f(x0[-1] + 3 / 8 * h, y0[-1] + 3 / 32 * k1 + 9 / 32 * k2)
			# k4 = h * f(x0[-1] + 12 / 13 * h, y0[-1] + 1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3)
			# k5 = h * f(x0[-1] + h, y0[-1] + 439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4)
			# k6 = h * f(x0[-1] + 1 / 2 * h, y0[-1] - 8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5)
			# y = y0[-1] + 16 / 135 * k1 + 6656 / 12825 * k3 + 28561 / 56430 * k4 - 9 / 50 * k5 + 2 / 55 * k6
			k0 = h * f(x0[-1], y0[-1])
			k1 = h * f(x0[-1] + 1 / 3 * h, y0[-1] + 1 / 3 * k0)
			k2 = h * f(x0[-1] + 1 / 3 * h, y0[-1] + 1 / 6 * k0 + 1 / 6 * k1)
			k3 = h * f(x0[-1] + 1 / 2 * h, y0[-1] + 1 / 8 * k0 + 3 / 8 * k2)
			k4 = h * f(x0[-1] + h, y0[-1] + 1 / 2 * k2 - 3 / 2 * k2 + 2 * k3)
			y = y0[-1] + (k0 + 4 * k3 + k4) / 6
			termination_status = 0
			h0 = h
		# print(y-y1)
		else:
			termination_status = 0
			# !!! Feedback to 'dt' in simulation required for correct work !!!
			tolerance = 1e-11
			check = False
			h0 = h
			# print(y0[-6:])
			while not check:
				y = (h0 * f(x0[-1], y0[-1]) - np.sum(y0[-6:] * a_predictor[:-1], axis=0, dtype='float64')) / a_predictor[6]
				e = tolerance1 + 1
				counter = 1
				while e >= tolerance1:
					y_prev = y
					y = (h0 * f(x0[-1] + h0, y) - np.sum(y0[-6:] * a_corrector[:-1], axis=0, dtype='float64')) / a_corrector[6]
					e = max(abs(y - y_prev))
					counter += 1
				print(min(y), end= '\t')
				k0 = h0 * f(x0[-1], y0[-1])
				k1 = h0 * f(x0[-1] + 1 / 3 * h0, y0[-1] + 1 / 3 * k0)
				k2 = h0 * f(x0[-1] + 1 / 3 * h0, y0[-1] + 1 / 6 * k0 + 1 / 6 * k1)
				k3 = h0 * f(x0[-1] + 1 / 2 * h0, y0[-1] + 1 / 8 * k0 + 3 / 8 * k2)
				k4 = h0 * f(x0[-1] + h0, y0[-1] + 1 / 2 * k2 - 3 / 2 * k2 + 2 * k3)
				# y = y0 + (k0 + 4 * k3 + k4) / 6
				R = (2 * k0 - 9 * k2 + 8 * k3 - k4) / 30
				# print(max(R))
				check1 = max(abs(R)) <= tolerance
				check2 = max(abs(R)) >= (tolerance / 30)
				'''Check if step status is changed by mass balances'''
				check3 = balance_check_inintegrator
				# print(max(abs(R)))
				# print('check1 - {}\tcheck2 - {}\tcheck3 - {}'.format(check1, check2, check3))
				if check1 and check2:
					break
				elif not check1:
					termination_status = 1
					return y, h0 * 0.75, termination_status
				# use h0 * 0.5
				elif not check2:
					if not check3:
						termination_status = 2
						return y, h0 * 1.5, termination_status
					else:
						break
		return y, h0, termination_status
			# цикл лекций лобанов аристова
			# химреактор конференция

		# print(e)
		# return y, h, 0


	else:
		print('ERROR! Specified integration method for function "{}" is not available. Specify valid'
			  'integration method or check spelling'.format(f.__name__))
		sys.exit()
