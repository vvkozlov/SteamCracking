"""
Header      : usermath.py
Created     : 08.01.2023
Author      : Vladimir Kozlov, kozlov.vlr@yandex.ru
Description : Contains some custom functions

References	:
			[1] - Add references for integration methods

"""
import sys


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


def integrate(method: str, f, x0: float, y0: float, h: float):
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
    else:
        print('ERROR! Specified integration method for function "{}" is not available. Specify valid'
              'integration method or check spelling'.format(f.__name__))
        sys.exit()
