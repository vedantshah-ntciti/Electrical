from ElectricalCalculation import ElectricalFormulae as ef
import math, cmath
import numpy as np

def solve():
    a= ef.linearequation(np.array([[2,-1],[-47,52]]),np.array([1,235]))
    return a[1]/47000

print(solve())