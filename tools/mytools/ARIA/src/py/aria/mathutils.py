"""

          ARIA -- Ambiguous Restraints for Iterative Assignment

                 A software for automated NOE assignment

                               Version 2.3


Copyright (C) Benjamin Bardiaux, Michael Habeck, Therese Malliavin,
              Wolfgang Rieping, and Michael Nilges

All rights reserved.


NO WARRANTY. This software package is provided 'as is' without warranty of
any kind, expressed or implied, including, but not limited to the implied
warranties of merchantability and fitness for a particular purpose or
a warranty of non-infringement.

Distribution of substantively modified versions of this module is
prohibited without the explicit permission of the copyright holders.

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""


from numpy import *

def average(x, n = None, exponent = 1., axis = 0):
    """
    Returns (n^{-1} sum_1^n x_i^exponent)^{1/exponent}.
    if 'n' is not None, it is used instead of len(x)
    sum is taken wrt to axis 'axis'
    """
    x = array(x, float)

    if n is None:
        n = shape(x)[axis]

    return (sum(power(x, exponent), axis) / n) ** (1. / exponent)

def _average(x):
    return sum(array(x), axis = 0) / len(x)

def variance(x, avg = None):
    if avg is None:
        avg = _average(x)

    return sum(power(array(x) - avg, 2), axis = 0) / (len(x) - 1.)

def standardDeviation(x, avg = None):
    return sqrt(variance(x, avg))

def confidenceInterval(x, p):
    """
    returns the smallest interval (start, end) and its size
    that covers at least a fraction of 'p' of the data-points
    given by x.
    """

    import math

    x = sort(array(x))

    ## set n to the next integer that is greater or equal
    ## to #data-points * p.

    n = len(x)
    m = n * p
    
    if math.ceil(n) > 0.:
        m = int(m) + 1
    else:
        m = int(m)

    smallest = x[-1] - x[0]
    start = x[0]
    stop = x[-1]

    for i in range(len(x)):

        if i + m >= len(x):
            interval = x[-1] - x[i] + x[m - n + i] - x[0]
            if interval < smallest:
                smallest = interval
                start = x[i]
                stop = x[m - n + i]

        else:
            interval = x[i + m] - x[i]
            if interval < smallest:
                smallest = interval
                start = x[i]
                stop = x[i + m]

    return smallest, start, stop
            
