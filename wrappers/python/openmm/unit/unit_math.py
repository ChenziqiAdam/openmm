"""
Module openmm.unit.math

Arithmetic methods on Quantities and Units

This is part of the OpenMM molecular simulation toolkit.
See https://openmm.org/development.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Christopher M. Bruns
Contributors: Peter Eastman

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import division, print_function, absolute_import

__author__ = "Christopher M. Bruns"
__version__ = "0.5"


import math
from .quantity import is_quantity
from .unit_definitions import *

try:
    from .. import _scientific_checkers as _scibench_checkers
except Exception:
    _scibench_checkers = None

####################
### TRIGONOMETRY ###
####################

def sin(angle):
    """
    Examples

    >>> sin(90*degrees)
    1.0
    """
    if is_quantity(angle):
        return math.sin(angle/radians)
    else:
        return math.sin(angle)

def sinh(angle):
    if is_quantity(angle):
        return math.sinh(angle/radians)
    else:
        return math.sinh(angle)

def cos(angle):
    """
    Examples

    >>> cos(180*degrees)
    -1.0
    """
    if is_quantity(angle):
        return math.cos(angle/radians)
    else:
        return math.cos(angle)

def cosh(angle):
    if is_quantity(angle):
        return math.cosh(angle/radians)
    else:
        return math.cosh(angle)

def tan(angle):
    if is_quantity(angle):
        return math.tan(angle/radians)
    else:
        return math.tan(angle)

def tanh(angle):
    if is_quantity(angle):
        return math.tanh(angle/radians)
    else:
        return math.tanh(angle)

def acos(x):
    """
    >>> acos(1.0)
    Quantity(value=0.0, unit=radian)
    >>> print(acos(1.0))
    0.0 rad
    """
    return math.acos(x) * radians

def acosh(x):
    return math.acosh(x) * radians

def asin(x):
    return math.asin(x) * radians

def asinh(x):
    return math.asinh(x) * radians

def atan(x):
    return math.atan(x) * radians

def atanh(x):
    return math.atanh(x) * radians

def atan2(x, y):
    return math.atan2(x, y) * radians

###################
### SQUARE ROOT ###
###################

def sqrt(val):
    """
    >>> sqrt(9.0)
    3.0
    >>> print(sqrt(meter*meter))
    meter
    >>> sqrt(9.0*meter*meter)
    Quantity(value=3.0, unit=meter)
    >>> sqrt(9.0*meter*meter*meter)
    Traceback (most recent call last):
    ...
    ArithmeticError: Exponents in Unit.sqrt() must be even.
    """
    try:
        result = val.sqrt()
    except AttributeError:
        return math.sqrt(val)
    if _scibench_checkers is not None and _scibench_checkers.enabled():
        try:
            if is_quantity(val) and isinstance(result.value_in_unit(result.unit), (int, float)):
                # OM-UNIT-005: squaring the result (via an independent
                # __mul__ path, not result.unit's own sqrt/conversion
                # machinery again) must reproduce the ORIGINAL quantity's
                # value, expressed in result.unit**2 -- not val._value
                # directly, since val may be expressed in a unit with a
                # nontrivial conversion factor to result.unit**2.
                got_val = result.value_in_unit(result.unit)
                squared = got_val * got_val
                orig_in_squared_unit = val.value_in_unit(result.unit ** 2)
                _scibench_checkers.check_sqrt_paths_agree(squared, orig_in_squared_unit, orig_in_squared_unit)
        except Exception:
            pass
    return result

###########
### SUM ###
###########

def sum(val):
    """
    >>> sum((1.0, 2.0))
    3.0
    >>> sum((2.0*meter, 3.0*meter))
    Quantity(value=5.0, unit=meter)
    >>> sum((2.0*meter, 30.0*centimeter))
    Quantity(value=2.3, unit=meter)
    """
    try:
        return val.sum()
    except AttributeError:
        pass
    if len(val) == 0:
        return 0
    result = val[0]
    for i in range(1, len(val)):
        result += val[i]
    return result

###################
### VECTOR MATH ###
###################

def dot(x, y):
    """
    >>> dot((2, 3)*meter, (4, 5)*meter)
    Quantity(value=23, unit=meter**2)
    """
    sum = x[0]*y[0]
    for i in range(1, len(x)):
        sum += x[i]*y[i]
    return sum

def norm(x):
    """
    >>> norm((3, 4)*meter)
    Quantity(value=5.0, unit=meter)
    """
    d = dot(x, x)
    result = sqrt(d)
    if _scibench_checkers is not None and _scibench_checkers.enabled():
        try:
            if is_quantity(x[0]):
                base_unit = x[0].unit
                dot_compat = is_quantity(d) and d.unit.is_compatible(base_unit ** 2)
                norm_compat = is_quantity(result) and result.unit.is_compatible(base_unit)
                _scibench_checkers.check_dot_norm_dimensions(dot_compat, norm_compat)
        except Exception:
            pass
    return result

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
