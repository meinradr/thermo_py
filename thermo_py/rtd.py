# TODO: make all numbers decimal.Decimal -> at the moment not possible since no way to use sympy lambdas with it
# TODO: to remove redundant code combine compare_rt and compare_t in one function
# TODO: create functions for printing values of single RTD similar to compare_rt
# TODO: write proper doc string for all methods, classes and module

from thermo_py.init import *
from thermo_py.init import _D0, _D1, _D2, _D3, _D100, _OHM


class RTD(object):
    _A = None
    _B = None
    _C = None
    _a = None
    _d = None
    _b = None
    _r0 = None
    _rs = None      # additional serial resistance to Rpt
    _min_t = -200
    _max_t = 850
    _min_rt = None
    _max_rt = None
    _tol_abs = None
    _tol_rel = None

    _f_rt = None            # name for function returning rt
    _f_t_pos = None         # name for function returning t above 0°C
    _f_zero = None          # name for function returning 0 if t and rt are correct
    _f_zero_d_t = None      # derivative of t for function _f_zero

    _zero = sympy.symbols('zero')
    _s_t, _s_rt, _s_r0, _s_A, _s_B, _s_C, _s_rs = sympy.symbols('t rt r0 A B C, rs')

    _eq = sympy.Eq(_zero,
                   _s_r0 * (_D1 + _s_A * _s_t + _s_B * _s_t ** _D2 + _s_C * (_s_t - _D100) * _s_t ** _D3) + _s_rs - _s_rt)

    _standard_n = None

    _standards = {
        # Coefficients from "Measuring Temperature with RTDs - A Tutorial" by National Instruments Corporation:
        # https://newton.ex.ac.uk/teaching/CDHW/Sensors/an046.pdf
    
        # TODO: Missing MIL-P-25726B and MIL-P-27723E referenced in Goodriche TAT report

        "EN 60751":{
            # DIN EN 60751 / IEC 751 standard
            'given' : 'ABC',
            'A' : '3.9083e-3',
            'B' :'-5.775e-7',
            'C' : '-4.183e-12'},

        "DIN43760":{
            # DIN 43760 standard.
            'given': 'ABC',
            'A' : '3.9080e-3',
            'B' : '-5.8019e-7',
            'C' : '-4.2735e-12'},

        "American":{
            # American standard.   ????ASTM E1137???
            'given': 'ABC',
            'A' : '3.9692e-3',
            'B' : '-5.8495e-7',
            'C' : '-4.2325e-12'},

        "ITS-90":{
            # ITS-90 standard.
            'given': 'ABC',
            'A' : '3.9848e-3',
            'B' : '-5.870e-7',
            'C' : '-4.0000e-12'},

        "SAMA":{
            'given': 'ABC',
            'A' : '3.97869E-3',
            'B' : '-5.86863E-7',
            'C' : '-4.16696E-12'},

        "IPTS-48":{
            # Goodrich TAT IPTS-48
            'given': 'adb',
            'a' : '0.003925',
            'd' : '1.45',
            'b' : '0.1'},

        "IPTS-68":{
            # Goodrich TAT IPTS-68
            'given': 'adb',
            'a' : '0.003925',
            'd' : '1.46',
            'b' : '0.1'}}

    _tol_class_n = None

    _tol_classes = {
        # TODO: EN 60751 Class F0.xx Number missing in comments

        # EN 60751 Class F0.10 (Class AA)  ±(0.1  + 0.0017 |t|)
        'EN 60751 AA': {
          'tol_abs': '0,1000',
          'tol_rel': '0.0017'},

        # EN 60751 Class F0.15 (Class A)   ±(0.15 + 0.002  |t|)
        'EN 60751 A': {
          'tol_abs': '0,1500',
          'tol_rel': '0.002'},

        # EN 60751 Class F0.x  (Class 1/3 B)   ±(0.3  + 0.005  |t|)
        'EN 60751 1/3 B': {
            'tol_abs': '0.1000',
            'tol_rel': '0.0050'},

        # EN 60751 Class F0.x  (Class 1/2 B)   ±(0.3  + 0.005  |t|)
        'EN 60751 1/2 B': {
            'tol_abs': '0.1500',
            'tol_rel': '0.0050'},

        # EN 60751 Class F0.3  (Class B)   ±(0.3  + 0.005  |t|)
        'EN 60751 B': {
          'tol_abs': '0.3000',
          'tol_rel': '0.0050'},

        # EN 60751 Class F0.x  (Class 2 B)   ±(0.3  + 0.005  |t|)
        'EN 60751 2 B': {
            'tol_abs': '0.6000',
            'tol_rel': '0.0050'},

        # EN 60751 Class F0.x  (Class 2 B)   ±(0.3  + 0.005  |t|)
        'EN 60751 5 B': {
            'tol_abs': '0.15000',
            'tol_rel': '0.0050'},

        # EN 60751 Class F0.6  (Class C)   ±(0.6  + 0.01   |t|)
        'EN 60751 C': {
          'tol_abs': '0.6000',
          'tol_rel': '0.0100'},

        # ASTM E1137 Grade A               ±(0.13 + 0.0017 |t|)
        'ASTM E1137 A': {
          'tol_abs': '0.1200',
          'tol_rel': '0.0017'},

        # ASTM E1137 Grade B               ±(0.25 + 0.0042 |t|)
        'ASTM E1137 B': {
          'tol_abs': '0.2500',
          'tol_rel': '0.0042'}}

    def __init__(self, r0=100, standard="EN 60751", tol_class="EN 60751 A", rs=0):
        self._create_functions()
        self.r0 = r0
        self.rs = rs
        self.standard = standard
        self.tolerance_class = tol_class

    @property
    def standard(self):
        return self._standard_n

    @standard.setter
    def standard(self, standard):
        if standard not in self._standards:
            raise ValueError(f'Unknown standard {standard}')

        s = self._standards[standard]

        if s['given'] == 'ABC':
            self.ABC = s['A'], s['B'], s['C']

        elif s['given'] == 'adb':
            self.adb = s['a'], s['d'], s['b']

        else:
            raise ValueError(f'No coefficients are specified in standard: {standard}')

        self._standard_n = standard

    @property
    def tolerance_class(self):
        return self._tol_class_n

    @tolerance_class.setter
    def tolerance_class(self, tol_class_n):
        if tol_class_n not in self._tol_classes:
            raise ValueError("Unknown tolerance class %s" % tol_class_n)

        self._tol_class_n = tol_class_n

    @property
    def name(self):
        return self._standard_n

    @name.setter
    def name(self, name):
        self._standard_n = name

    @property
    def ABC(self):
        return self._A, self._B, self._C

    @ABC.setter
    def ABC(self, ABC):
        self._A = Decimal(ABC[0])
        self._B = Decimal(ABC[1])
        self._C = Decimal(ABC[2])

        self._a, self._d, self._b = self.ABC_to_adb(self._A, self._B, self._C)

        # print(f'{self.A}   {self.B}   {self.C}    {self.a}   {self.d}   {self.b}')

        self._reset_min_max_rt()

        self._standard_n = None

    @property
    def adb(self):
        return self._a, self._d, self._b

    @adb.setter
    def adb(self, adb):

        self._a = Decimal(adb[0])
        self._d = Decimal(adb[1])
        self._b = Decimal(adb[2])

        self._A, self._B, self._C = self.adb_to_ABC(self._a, self._d, self._b)

        #print(f'{self.A}   {self.B}   {self.C}    {self.a}   {self.d}   {self.b}')

        self._reset_min_max_rt()

        self._standard_n = None

    @property
    def r0(self):
        return self._r0

    @r0.setter
    def r0(self, r0):
        self._r0 = r0

        self._reset_min_max_rt()

    @property
    def rs(self):
        return self._rs

    @rs.setter
    def rs(self, rs):
        self._rs = rs

        self._reset_min_max_rt()

    def _reset_min_max_rt(self):
        if self._r0 and self._A and self._B and self._C:
            self._min_rt = self.calc_rt(self._min_t) - 1
            self._max_rt = self.calc_rt(self._max_t) + 1

    def get_tol(self, t, tol_class_n = None):

        if tol_class_n:
            if tol_class_n not in self._tol_classes:
                raise ValueError("Unknown tolerance class %s" % tol_class_n)
        else:
            if self._tol_class_n:
                tol_class_n = self._tol_class_n
            else:
                raise ValueError("Tolerance class has not been defined")

        tol_class = self._tol_classes[tol_class_n]

        tol_abs = Decimal(tol_class['tol_abs'])
        tol_rel = Decimal(tol_class['tol_abs'])

        tol = tol_abs + abs(t) * tol_rel

        return tol

    def _create_functions(self):
        """ creating functions for calculating rt and t and some helper functions"""

        eq = self._eq

        t = self._s_t
        rt = self._s_rt
        r0 = self._s_r0
        A = self._s_A
        B = self._s_B
        C = self._s_C
        zero = self._zero
        rs = self._s_rs

        # ---- creating _f_rt -------------------------------------------------
        # replacing _zero symbol with value 0
        eq_rt = eq.subs(zero, 0)
        # solve for rt
        solutions_rt = sympy.solve(eq_rt, rt)
        # take first solution
        s_rt = solutions_rt[0]
        # create function
        self._f_rt = sympy.lambdify([t, r0, A, B, C, rs], s_rt)

        # ---- creating _f_t_pos ----------------------------------------------
        # replacing _zero and C symbol with value 0
        eq_t_pos = eq.subs([(zero, 0), (C, 0)])
        # solve for t
        solutions_t_pos = sympy.solve(eq_t_pos, t)
        # take first solution
        s_t_pos = solutions_t_pos[0]
        # create function
        self._f_t_pos = sympy.lambdify([rt, r0, A, B, rs], s_t_pos)

        # ---- creating _f_zero ----------------------------------------------
        # solve for zero
        solutions_zero = sympy.solve(eq, zero)
        # take first solution
        s_zero = solutions_zero[0]
        # create function
        self._f_zero = sympy.lambdify([t, rt, r0, A, B, C, rs], s_zero)

        # ---- creating _f_zero_d_t ----------------------------------------------
        # differentiate solution of _f_zero
        s_zero_d_t = sympy.diff(s_zero, t)
        # create function
        self._f_zero_d_t = sympy.lambdify([t, r0, A, B, C], s_zero_d_t)

    def calc_rt(self, t):
        """ calculating rt from t"""

        # cast input to decimal
        if isinstance(t, str):
            t = Decimal(t)

        # check if input is within the specified temperature range
        if not self._min_t <= t <= self._max_t:
            raise ValueError(f'temperature {t}°C is out of range'
                             f'({self._min_t}°C to {self._max_t}°C) for an RTD')

        # get some variables
        r0 = self._r0
        A = self._A
        B = self._B
        rs = self._rs

        # set C to = if input is > 0
        if t < 0:
            C = self._C
        else:
            C = _D0

        # convert to float since sympy lambdas do not support Decimal
        t = float(t)
        r0 = float(r0)
        A = float(A)
        B = float(B)
        C = float(C)
        rs = float(rs)


        # calculate rt
        rt = self._f_rt(t, r0, A, B, C, rs)

        return rt

    def calc_t(self, rt, tolerance=1e-6):

        rt = Decimal(rt)

        # check if rt is in range of specified RTD
        if not self._min_rt <= rt <= self._max_rt:
            raise ValueError(f'resistance (rt) {rt} {_OHM} is out of range'
                             f'({self._min_rt} {_OHM} to {self._max_rt} {_OHM}) for the specified RTD')

        # get some variables
        r0 = self._r0
        A = self._A
        B = self._B
        C = self._C
        rs = self._rs

        # convert to float since sympy lambdas do not support Decimal
        rt = float(rt)
        r0 = float(r0)
        A = float(A)
        B = float(B)
        C = float(C)
        rs = float(rs)

        t = self._f_t_pos(rt, r0, A, B, rs)

        # calculate t with the formula including C for negative numbers
        if rt < r0:
            # calculate t with Newton method
            f = lambda t_: self._f_zero(t_, rt, r0, A, B, C, rs)
            Df = lambda t_: self._f_zero_d_t(t_, r0, A, B, C)
            t = self.newton_method(f, Df, t, tolerance, 10)

        return t

    @staticmethod
    def ABC_to_adb(A, B, C=0):
        """ Calculates Alpha, Delta and Beta form the A, B and C coenficens
            source: 1399021951CalVan.pdf"""
        a = A - (-100 * B)
        d = (-1 * B * (100 ** 2)) / a
        b = (-1 * C * (100 ** 4)) / a

        return a, d, b

    @staticmethod
    def adb_to_ABC(a, d, b=0):
        """ Calculates the A, B and C coenficens form Alpha, Delta and Beta
            source: 1399021951CalVan.pdf"""

        A = a + ((a * d) / 100)
        B = -1 * ((a * d) / (100 ** 2))
        C = -1 * ((a * b) / (100 ** 4))

        return A, B, C

    @staticmethod
    def calc_alpha(R0, R100):
        # TODO: to be double checked

        # RTD alpha coefficient
        # R0 resistance at 0C
        # R100 resistance at 100C
        # http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
        return (R100 - R0) / (100 * R0)

    @staticmethod
    def calc_delta(r0, Rth, Th, alpha):
        # TODO: to be double checked

        # RTD delta coefficient
        # R0 resistance at 0C
        # Th highest temperature in the calibration range
        # Rth resistance of the sensor at the highest temperature
        # John G. Webster and Halit Eren, 2014
        # Measurement, Instrumentation, and Sensors Handbook, Second Edition
        # Spatial, Mechanical, Thermal, and Radiation Measurement
        # CRC Press

        return (Th - (Rth - r0) / (r0 * alpha)) / ((Th / 100 - 1) * (Th / 100))

    @staticmethod
    def calc_beta(R0, Rtl, Tl, alpha, delta):
        # TODO: to be double checked

        # RTD beta coefficient
        # R0 resistance at 0C
        # Tl lowest temperature in the calibration range
        # Rtl resistance of the sensor at the lowest temperature
        # John G. Webster and Halit Eren, 2014
        # Measurement, Instrumentation, and Sensors Handbook, Second Edition
        # Spatial, Mechanical, Thermal, and Radiation Measurement
        # CRC Press

        return (Tl - ((Rtl - R0) / (R0 * alpha) + delta * (Tl / 100 - 1) * (Tl / 100))) / ((Tl / 100 - 1) * (Tl / 100) ** 3)

    @staticmethod
    def newton_method(f, Df, x0, epsilon, max_iter):
        '''Approximate solution of f(x)=0 by Newton's method.

        Parameters
        ----------
        f : function
            Function for which we are searching for a solution f(x)=0.
        Df : function
            Derivative of f(x).
        x0 : number
            Initial guess for a solution f(x)=0.
        epsilon : number
            Stopping criteria is abs(f(x)) < epsilon.
        max_iter : integer
            Maximum number of iterations of Newton's method.

        Returns
        -------
        xn : number
            Implement Newton's method: compute the linear approximation
            of f(x) at xn and find x intercept by the formula
                x = xn - f(xn)/Df(xn)
            Continue until abs(f(xn)) < epsilon and return xn.
            If Df(xn) == 0, return None. If the number of iterations
            exceeds max_iter, then return None.

        Examples
        --------
        >> f = lambda x: x**2 - x - 1
        >> Df = lambda x: 2*x - 1
        >> newton(f,Df,1,1e-8,10)
        Found solution after 5 iterations.
        1.618033988749989

        source: https://www.math.ubc.ca/~pwalls/math-python/roots-optimization/newton/
        '''

        xn = x0
        for n in range(0, max_iter):
            fxn = f(xn)
            if abs(fxn) < epsilon:
                #print('Found solution after', n, 'iterations.')
                return xn
            Dfxn = Df(xn)
            if Dfxn == 0:
                # print('Zero derivative. No solution found.')
                raise NoSolutionException(f'Newton Method: no solution found for x0: {x0} - Zero derivative')
            xn = xn - fxn / Dfxn
        # print('Exceeded maximum iterations. No solution found.')
        raise NoSolutionException(f'Newton Method: no solution found for x0: {x0} - Exceeded maximum iterations {max_iter}')

    @staticmethod
    def compare_rt(rtd_1, rtd_2, t_values, csv='', precision=2):
        p = precision
        result = ''
        result += '\n'

        result += f'{f"%{p + 12}s" % "name"}' + csv
        result += f'{f"%{p + 12}s" % rtd_1.standard}' + csv
        result += f'{f"%{p + 12}s" % rtd_2.standard}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "r0"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_1.r0)}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_2.r0)}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "A"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_1.ABC[0])}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_2.ABC[0])}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "B"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_1.ABC[1])}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_2.ABC[1])}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "C"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_1.ABC[2])}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_2.ABC[2])}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "Rs"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_1.rs)}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_2.rs)}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "t [°C]"}' + csv
        result += f'{f"%{p + 12}s" % ("Rt 1 [" + _OHM + "]")}' + csv
        result += f'{f"%{p + 12}s" % ("Rt 2 [" + _OHM + "]")}' + csv
        result += f'{f"%{p + 12}s" % ("dev [" + _OHM + "]")}' + csv
        result += '\n'

        for t in t_values:

            rt_1 = rtd_1.calc_rt(t)
            rt_2 = rtd_2.calc_rt(t)
            err = rt_1 - rt_2

            result += f'{f"%{p + 12}s" % (f"%.{p}f" % t)}' + csv
            result += f'{f"%{p + 12}s" % (f"%.{p}f" % rt_1)}' + csv
            result += f'{f"%{p + 12}s" % (f"%.{p}f" % rt_2)}' + csv
            result += f'{f"%{p + 12}s" % (f"%.{p}f" % err)}' + csv
            result += '\n'

        return result

    @staticmethod
    def compare_t(rtd_1, rtd_2, rt_values, csv='', precision=2):
        p = precision
        result = ''

        result += f'{f"%{p + 12}s" % "name"}' + csv
        result += f'{f"%{p + 12}s" % rtd_1.standard}' + csv
        result += f'{f"%{p + 12}s" % rtd_2.standard}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "r0"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_1.r0)}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_2.r0)}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "A"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_1.ABC[0])}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_1.ABC[0])}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "B"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_1.ABC[1])}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_2.ABC[1])}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "C"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_1.ABC[2])}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 2}E" % rtd_2.ABC[2])}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % "Rs"}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_1.rs)}' + csv
        result += f'{f"%{p + 12}s" % (f"%.{p + 1}f" % rtd_2.rs)}' + csv
        result += '\n'

        result += f'{f"%{p + 12}s" % ("Rt [" + _OHM + "]")}' + csv
        result += f'{f"%{p + 12}s" % "t 1 [°C]"}' + csv
        result += f'{f"%{p + 12}s" % "t 2 [°C]"}' + csv
        result += f'{f"%{p + 12}s" % "dev [°C]"}' + csv
        result += '\n'

        for rt in rt_values:

            t_1 = rtd_1.calc_t(rt)
            t_2 = rtd_2.calc_t(rt)
            err = t_1 - t_2

            result += f'{f"%{p + 12}s" % (f"%.{p}f" % rt)}' + csv
            result += f'{f"%{p + 12}s" % (f"%.{p}f" % t_1)}' + csv
            result += f'{f"%{p + 12}s" % (f"%.{p}f" % t_2)}' + csv
            result += f'{f"%{p + 12}s" % (f"%.{p}f" % err)}' + csv
            result += '\n'

        return result
