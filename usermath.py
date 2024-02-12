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

import coreobjects


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
		return y0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), h, 0
	elif method == 'rungekuttafelberg5th':
		k1 = h * f(x0, y0)
		k2 = h * f(x0 + 1 / 4 * h, y0 + 1 / 4 * k1)
		k3 = h * f(x0 + 3 / 8 * h, y0 + 3 / 32 * k1 + 9 / 32 * k2)
		k4 = h * f(x0 + 12 / 13 * h, y0 + 1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3)
		k5 = h * f(x0 + h, y0 + 439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4)
		k6 = h * f(x0 + 1 / 2 * h, y0 - 8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5)
		return y0 + 16 / 135 * k1 + 6656 / 12825 * k3 + 28561 / 56430 * k4 - 9 / 50 * k5 + 2 / 55 * k6, h, 0
	elif method == 'rungekuttamerson_fixed':
		h0 = h
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
		return y, h, 0, R
	elif method == 'rungekuttamerson_adaptive':
		'''termination_status stores the reason of integration termination:
			0 - good accuracy
			1 - error is too large, reduce dh
			2 - error is too small, increase dh'''
		termination_status = 0
		# !!! Feedback to 'dt' in simulation required for correct work !!!
		truncation_tolerance = 1e-11
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
			check1 = R <= truncation_tolerance
			check2 = R >= (truncation_tolerance / 30)
			'''Check if step status is changed by mass balances'''
			check3 = balance_check_inintegrator
			# print(max(abs(R)))
			# print('check1 - {}\tcheck2 - {}'.format(check1, check2))
			if check1 and check2:
				break
			elif not check1:
				termination_status = 1
				return y, h0 * 0.75, termination_status, R
			# use h0 * 0.5
			elif not check2:
				if not check3:
					termination_status = 2
					return y, h0 * 1.5, termination_status, R
				else:
					break
		return y, h0, termination_status, R
		"""
		function y = runge(t0, y0, dt, h, tol)
		
		flag = 1;
		t = t0;
		y = y0;
		while (t < t0 + dt)
		    k1 = h * f(t, y);
		    k2 = h * f(t + 1/4*h,   y + 1/4*k1);
		    k3 = h * f(t + 3/8*h,   y + 3/32*k1 +       9/32*k2);
		    k4 = h * f(t + 12/13*h, y + 1932/2197*k1 -  7200/2197*k2 +  7296/2197*k3);
		    k5 = h * f(t + h,       y + 439/216*k1 -    8*k2 +          3680/513*k3 - 845/4104*k4);
		    k6 = h * f(t + h/2,     y - 8/27*k1 +   2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5);
		    if (flag)
		        z = y + 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6;
		    end
		    y = y + 25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - k5/5;
		    if (flag && abs(y - z) > tol)
		        t = t0;
		        y = y0;
		        h = h * (tol * h/2/abs(z-y))^0.25;
		    else
		        t = t + h;
		        flag = 0;
		    end
		end
		"""
	elif method == 'gear':
		corrector_tolerance = 1e-14
		h0 = h
		if len(x0) < 60:
			termination_status = 0
			# !!! Feedback to 'dt' in simulation required for correct work !!!
			truncation_tolerance = 1e-11
			k0 = h * f(x0[-1], y0[-1])
			k1 = h * f(x0[-1] + 1 / 3 * h, y0[-1] + 1 / 3 * k0)
			k2 = h * f(x0[-1] + 1 / 3 * h, y0[-1] + 1 / 6 * k0 + 1 / 6 * k1)
			k3 = h * f(x0[-1] + 1 / 2 * h, y0[-1] + 1 / 8 * k0 + 3 / 8 * k2)
			k4 = h * f(x0[-1] + h, y0[-1] + 1 / 2 * k2 - 3 / 2 * k2 + 2 * k3)
			y = y0[-1] + (k0 + 4 * k3 + k4) / 6
			"""
			check = False
			while not check:
				### FIRST STEPS (in merson) ARE NOT WORKING PROPERLY!
				
				R = (2 * k0 - 9 * k2 + 8 * k3 - k4) / 30
				try:
					R = max(abs(R))
				except:
					R = abs(R)
				# print(R)
				check1 = R <= truncation_tolerance
				check2 = R >= (truncation_tolerance / 30)
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
			"""
			# k1 = h * f(x0, y0[-1])
			# k2 = h * f(x0[-1] + 1 / 4 * h, y0[-1] + 1 / 4 * k1)
			# k3 = h * f(x0[-1] + 3 / 8 * h, y0[-1] + 3 / 32 * k1 + 9 / 32 * k2)
			# k4 = h * f(x0[-1] + 12 / 13 * h, y0[-1] + 1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3)
			# k5 = h * f(x0[-1] + h, y0[-1] + 439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4)
			# k6 = h * f(x0[-1] + 1 / 2 * h, y0[-1] - 8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5)
			# y = y0[-1] + 16 / 135 * k1 + 6656 / 12825 * k3 + 28561 / 56430 * k4 - 9 / 50 * k5 + 2 / 55 * k6
			"""
			k0 = h * f(x0[-1], y0[-1])
			k1 = h * f(x0[-1] + 1 / 3 * h, y0[-1] + 1 / 3 * k0)
			k2 = h * f(x0[-1] + 1 / 3 * h, y0[-1] + 1 / 6 * k0 + 1 / 6 * k1)
			k3 = h * f(x0[-1] + 1 / 2 * h, y0[-1] + 1 / 8 * k0 + 3 / 8 * k2)
			k4 = h * f(x0[-1] + h, y0[-1] + 1 / 2 * k2 - 3 / 2 * k2 + 2 * k3)
			y = y0[-1] + (k0 + 4 * k3 + k4) / 6
			termination_status = 0
			h0 = h
			"""
		# print(y-y1)
		else:
			termination_status = 0
			# !!! Feedback to 'dt' in simulation required for correct work !!!
			truncation_tolerance = 1e-11
			check = False
			h0 = h
			# print(y0[-6:])
			while not check:
				y = (h0 * f(x0[-1], y0[-1]) - np.sum(y0[-6:] * a_predictor[:-1], axis=0, dtype='float64')) / a_predictor[6]
				e = corrector_tolerance + 1
				counter = 1
				# print('\t\te0', e)
				while e >= corrector_tolerance:
					y_prev = y
					y = (h0 * f(x0[-1] + h0, y) - np.sum(y0[-6:] * a_corrector[:-1], axis=0, dtype='float64')) / a_corrector[6]
					e = max(abs(y - y_prev))
					# print('\te =', e)
					counter += 1
				# print(min(y), end= '\t')
				k0 = h0 * f(x0[-1], y0[-1])
				k1 = h0 * f(x0[-1] + 1 / 3 * h0, y0[-1] + 1 / 3 * k0)
				k2 = h0 * f(x0[-1] + 1 / 3 * h0, y0[-1] + 1 / 6 * k0 + 1 / 6 * k1)
				k3 = h0 * f(x0[-1] + 1 / 2 * h0, y0[-1] + 1 / 8 * k0 + 3 / 8 * k2)
				k4 = h0 * f(x0[-1] + h0, y0[-1] + 1 / 2 * k2 - 3 / 2 * k2 + 2 * k3)
				# y_runge = y0 + (k0 + 4 * k3 + k4) / 6
				R = (2 * k0 - 9 * k2 + 8 * k3 - k4) / 30
				# print('R = ', max(abs(R)))
				# print(max(R))
				check1 = max(abs(R)) <= truncation_tolerance
				check2 = max(abs(R)) >= (truncation_tolerance / 30)
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


def mbal_check(method: str, inlet_stream: coreobjects.Stream, outlet_stream: coreobjects.Stream, tolerance: float)\
		-> tuple[coreobjects.Stream, int]:
	'''Function descrription'''
	iterations = 0
	'''Initial balance check to be able to bypass function if everything is ok'''
	check = abs(inlet_stream.FLMASS - outlet_stream.FLMASS) <= tolerance
	if method == 'bisection':
		'''Initializing correctors'''
		corrector_l = 0.5
		corrector_r = 2
		'''Additional coorector variables to store previous corrector value to eliminate simulation call if
		                border did not change'''
		corrector_l_prev = []
		corrector_r_prev = []
		outlet_l = outlet_m = outlet_r = outlet_stream
		while not check:
			iterations += 1
			'''Calculate function value on the left border only if border has changed'''
			if corrector_l != corrector_l_prev or corrector_l_prev == []:
				outlet_l = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR,
											  outlet_stream.FLMOL * corrector_l, outlet_stream.P, outlet_stream.T, 'IG')
				corrector_l_prev = corrector_l
			else:
				pass
			'''Calculate new bisection point and function there'''
			corrector_m = (corrector_l + corrector_r) / 2
			outlet_m = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR,
										  outlet_stream.FLMOL * corrector_m, outlet_stream.P, outlet_stream.T, 'IG')
			'''Calculate function value on the right border only if border has changed'''
			if corrector_r != corrector_r_prev or corrector_r_prev == []:
				outlet_r = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR,
											  outlet_stream.FLMOL * corrector_r, outlet_stream.P, outlet_stream.T, 'IG')
			else:
				pass
			'''Calculate errors and s shift borders if required'''
			check_l = outlet_l.FLMASS - inlet_stream.FLMASS
			check_m = outlet_m.FLMASS - inlet_stream.FLMASS
			check_r = outlet_r.FLMASS - inlet_stream.FLMASS
			if abs(check_m) <= tolerance:
				check = True
			else:
				if check_m * check_l < 0:
					corrector_r = corrector_m
				elif check_m * check_r <0:
					corrector_l = corrector_m
				else:
					print('\t\tERROR! Check variable is equal to zero when trying to converge mass balance')
			'''Returning converged outlet stream'''
			# outlet_stream = outlet_m
		if not check:
			print('\t\tERROR! Mass balance has not been converged for some reason')
	elif method == 'bisection_gr':
		'''Initializing correctors'''
		corrector_l = 0.5
		corrector_r = 2
		golden_ratio = (5**0.5 + 1) / 2
		'''Additional coorector variables to store previous corrector value to eliminate simulation call if
						border did not change'''
		corrector_l_prev = []
		corrector_r_prev = []
		outlet_l = outlet_m = outlet_r = outlet_stream
		while not check:
			iterations += 1
			'''Calculate function value on the left border only if border has changed'''
			if corrector_l != corrector_l_prev or corrector_l_prev != []:
				outlet_l = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR,
											  outlet_stream.FLMASS * corrector_l, outlet_stream.P, outlet_stream.T,
											  'IG')
				corrector_l_prev = corrector_l
			else:
				pass
			'''Calculate new bisection point using golden ratio and function there'''
			corrector_m = corrector_l + (corrector_r - corrector_l) / golden_ratio
			outlet_m = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR,
										  outlet_stream.FLMASS * corrector_m, outlet_stream.P, outlet_stream.T, 'IG')
			'''Calculate function value on the right border only if border has changed'''
			if corrector_r != corrector_r_prev or corrector_r_prev != []:
				outlet_r = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR,
											  outlet_stream.FLMASS * corrector_r, outlet_stream.P, outlet_stream.T,
											  'IG')
			else:
				pass
			'''Calculate errors and s shift borders if required'''
			check_l = outlet_l.FLMASS - inlet_stream.FLMASS
			check_m = outlet_m.FLMASS - inlet_stream.FLMASS
			check_r = outlet_r.FLMASS - inlet_stream.FLMASS
			if abs(check_m) <= tolerance:
				check = True
			else:
				if check_m * check_l < 0:
					corrector_r = corrector_m
				elif check_m * check_r < 0:
					corrector_l = corrector_m
				else:
					print('\t\tERROR! Check variable is equal to zero when trying to converge mass balance')
			'''Returning converged outlet stream'''
			if not check:
				print('\t\tERROR! Mass balance has not been converged for some reason')
	elif method == 'secant':
		'''Initializing correctors'''
		x0 = 0.5
		x1 = 1.5
		'''Conditions for loop. Definitions of "tolerance" are a little bit messed up here'''
		check = abs(outlet_stream.FLMASS - inlet_stream.FLMASS) <= tolerance
		if check:
			iterations = 0
			pass
		else:
			error = abs(x0 - x1)
			'''Function under investigation is f(x) = outlet_stream.FLMASS - inlet_stream.FLMASS'''
			while not check:
				iterations += 1
				f_x0 = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR, outlet_stream.FLMASS * x0,
										  outlet_stream.P, outlet_stream.T, 'IG').FLMASS - inlet_stream.FLMASS
				f_x1 = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR, outlet_stream.FLMASS * x1,
										  outlet_stream.P, outlet_stream.T, 'IG').FLMASS - inlet_stream.FLMASS
				x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
				error = abs(x2 - x1)
				check = error < tolerance
				x0, x1 = x1, x2
			outlet_stream = coreobjects.Stream(outlet_stream.compset, outlet_stream.COMPMOLFR, outlet_stream.FLMASS * x2,
										  outlet_stream.P, outlet_stream.T, 'IG')
	else:
		check = False
		print('\t\tERROR! Selected method for MBAL convergence is not supported')
		sys.exit()
	return outlet_stream, iterations
