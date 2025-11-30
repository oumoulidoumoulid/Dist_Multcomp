from thermo.chemical import Chemical
import numpy as np

class Compound:
    """
    Represente un compose chimique avec ses proprietes
    """

    def __init__(self, name):
        self.name = name
        self.chem = Chemical(name)

        # Proprietes critiques
        self.Tc = self.chem.Tc  # Temperature critique (K)
        self.Pc = self.chem.Pc  # Pression critique (Pa)
        self.omega = self.chem.omega  # Facteur acentrique

        # Proprietes normales
        self.Tb = self.chem.Tb  # Temperature d'ebullition (K)
        self.MW = self.chem.MW  # Masse molaire (g/mol)

    def vapor_pressure(self, T):
        """ Pression de vapeur saturante a T (Pa)"""
        self.chem.T = T
        return self.chem.Psat

    def K_value(self, T, P):
        """ Coefficient de partage K = Psat / P"""
        Psat = self.vapor_pressure(T)
        return Psat / P if Psat else 0.0

    def enthalpy_liquid(self, T):
        """ Enthalpie du liquide a T (J/mol)"""
        self.chem.T = T
        self.chem.calculate(T)
        # Try different attributes for liquid enthalpy
        try:
            return self.chem.Hl
        except:
            try:
                return self.chem.H
            except:
                # Fallback: estimate using heat of vaporization
                # H_liquid ≈ -Hvap (rough approximation)
                try:
                    return -self.chem.Hvap if self.chem.Hvap else 0.0
                except:
                    return 0.0

    def enthalpy_vapor(self, T):
        """ Enthalpie de la vapeur a T (J/mol )"""
        self.chem.T = T
        self.chem.calculate(T)
        # Vapor enthalpy = Liquid enthalpy + Heat of vaporization
        try:
            H_l = self.enthalpy_liquid(T)
            H_vap = self.chem.Hvap if self.chem.Hvap else 0.0
            return H_l + H_vap
        except:
            try:
                return self.chem.H
            except:
                return 0.0
