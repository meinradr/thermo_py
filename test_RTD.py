import unittest
import json
from thermo_py import *


class TestRTD(unittest.TestCase):

    def test_omega_pt100_rt(self):
        tolerance = 0.015   # Ohm
        with open('thermo_py/OMEGA_PT100_Values.json') as json_file:
            data = json.load(json_file)['omega_pt100']
        r0 = data['r0']
        standard = data['standard']

        rtd = RTD(r0=r0, standard=standard)

        test_values = data['values']

        for pair in test_values:
            rt = rtd.calc_rt(pair['t'])
            err = abs(rt - pair['r'])
            try:
                self.assertTrue(err < tolerance)
            except AssertionError:
                raise AssertionError(f'err: {err} is larger than tolerance {tolerance} for t: {pair["t"]}')

    def test_omega_pt100_t(self):
        tolerance = 0.04    # °C
        with open('thermo_py/OMEGA_PT100_Values.json') as json_file:
            data = json.load(json_file)['omega_pt100']
        r0 = data['r0']
        standard = data['standard']

        rtd = RTD(r0=r0, standard=standard)

        test_values = data['values']

        for pair in test_values:
            t = rtd.calc_t(pair['r'])
            err = abs(t - pair['t'])
            self.assertTrue(err < tolerance)

    def test_calc_rt_vs_calc_t(self):
        tolerance = 1e-5
        factor = 10

        rtd = RTD()

        min = int(rtd._min_t) * factor
        max = int(rtd._max_t) * factor

        for i in range(min, max, 1):
            t = i / factor
            rt = rtd.calc_rt(t)
            t_ = rtd.calc_t(rt)
            err = abs(t - t_)
            self.assertTrue(err < tolerance)

    def test_ABD_adb_conversion(self):
        tolerance = 1e-29
        rtd = RTD()

        A, B, C = rtd.ABC
        a_, d_, b_ = rtd.ABC_to_adb(A, B, C)
        A_, B_, C_ = rtd.adb_to_ABC(a_, d_, b_)
        a, d, b = rtd.adb

        self.assertTrue(abs(A / A_ - 1) < tolerance)
        self.assertTrue(abs(B / B_ - 1) < tolerance)
        self.assertTrue(abs(C / C_ - 1) < tolerance)

        self.assertTrue(abs(a / a_ - 1) < tolerance)
        self.assertTrue(abs(d / d_ - 1) < tolerance)
        self.assertTrue(abs(b / b_ - 1) < tolerance)

    def test_MIL_formula(self):
        tolerance = 1e-12

        r0 = 100
        a = 0.003925
        d = 1.45
        b = 0.1

        rtd = RTD()
        rtd.r0 = r0
        rtd.adb = (a, d, b)

        def MIL_formula(t_, r0_, a_, d_, b_):
            if t_ >= 0:
                b_ = 0

            f1 = d_ * ((t_ / 100) - 1) * (t_ / 100)
            f2 = b_ * ((t_ / 100) - 1) * ((t_ / 100) ** 3)
            f3 = 1 + (a_ * (t_ - f1 - f2))
            return f3 * r0_

        for i in range(-200, 851, 1):
            mil_rt = MIL_formula(i, r0, a, d, b)
            rt = rtd.calc_rt(i)
            self.assertTrue(abs(mil_rt - rt) < tolerance)


if __name__ == '__main__':
    unittest.main()
