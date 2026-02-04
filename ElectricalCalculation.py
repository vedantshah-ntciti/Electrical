import numpy as np
import cmath
import math

class ElectricalFormulae:
    @staticmethod
    def series(*args:complex|float) -> float|complex:
        Net = sum(x for x in args)
        return Net
    
    @staticmethod
    def parallel(*args:complex|float) -> float|complex:
        if any(args == 0 for i in args):
            return 0
        Net = sum(1/x for x in args)
        return 1/Net
    
    @staticmethod
    def linearequation(A:np.array,B:np.array) -> np.array:
        if np.linalg.det(A)!=0:
            v=(np.linalg.solve(A,B))
            return v
        else:
            return np.nan
    
    @staticmethod
    def star_to_delta(ra, rb, rc) -> tuple[float,float,float]:
        """
        Converts Star (Wye) network resistances to Delta equivalent.
        """
        # Transformation formulas
        rab = ra + rb + (ra * rb / rc)
        rbc = ra + rb + (ra * rb / ra)
        rca = rb + ra + (rb * ra / ra)
        
        return rab,rbc,rca
    
    @staticmethod
    def delta_to_star(rab, rbc, rca) -> tuple[float,float,float]:
        """
        Converts Delta network resistances to Star (Wye) equivalent.
        """
        denom = rab + rbc + rca
        if denom == 0:
            return "Error: Sum of resistances must be greater than zero."
        
        # Transformation formulas
        ra = (rab * rca) / denom
        rb = (rab * rbc) / denom
        rc = (rbc * rca) / denom

        return ra,rb,rc
    
    @staticmethod
    def MillmansLaw(V:np.array,R:np.array) -> tuple[float,float]:
        if np.any(R==0):
            raise ValueError("Internal resistance cannot be zero")
        
        R_inv=1/R
        R_eff=1/np.sum(R_inv)
        E_eff=np.sum(V*R_inv)*R_eff

        return E_eff, R_eff

    @staticmethod
    def quadratic(a:complex|float,b:complex|float,c:complex|float) -> tuple[complex,complex]:
        if a == 0:
            return cmath.nan,cmath.nan
        
        d=cmath.sqrt(b*b - 4*a*c)
        
        sol1 = (-b+d)/(2*a)
        sol2 = (-b-d)/(2*a)
        return sol1,sol2