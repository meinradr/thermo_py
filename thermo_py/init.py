
# -------- importing modules for whole project -------------------------------------------------------

from decimal import Decimal
import sympy
from decimal import Decimal, getcontext


# -------- setting variables -------------------------------------------------------

# set Decimal precision
getcontext().prec = 30


# -------- defining constants -------------------------------------------------------------------

# defining Decimal constants
_D0 = Decimal('0')
_D1 = Decimal('1')
_D2 = Decimal('2')
_D3 = Decimal('3')
_D100 = Decimal('100')

# Constant for the greek Omega (Ohm) character
_OHM = '\u03A9'


# -------- defining Exceptions -------------------------------------------------------------------

class NoSolutionException(Exception):
    """ Exception if there is no solution in calculating t from rt"""
    pass
