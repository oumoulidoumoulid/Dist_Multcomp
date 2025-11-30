import numpy as np
from scipy.optimize import fsolve

class ThermodynamicPackage:
    """
    Calcule les proprietes thermodynamiques du melange
    """

    def __init__(self, compounds):
        """
        Parameters :
        -----------
        compounds : list of Compound objects
        """
        self.compounds = compounds
        self.n_comp = len(compounds)

    def K_values(self, T, P, x=None):
        """
        Calcule tous les K-values a T, P

        Returns :
        --------
        K : array de taille n_comp
        """
        K = np.zeros(self.n_comp)
        for i, comp in enumerate(self.compounds):
            K[i] = comp.K_value(T, P)
        return K

    def bubble_temperature(self, P, x, T_guess=None):
        """
        Calcule la temperature de bulle

        Returns :
        --------
        T_bubble : float (K)
        """
        if T_guess is None:
            # Estimation initiale
            T_guess = np.sum([x[i] * self.compounds[i].Tb for i in range(self.n_comp)])

        def equation(T):
            # fsolve passes T as an array, take the first element
            T_val = T[0] if isinstance(T, np.ndarray) else T
            K = self.K_values(T_val, P)
            return np.sum(K * x) - 1.0

        T_bubble = fsolve(equation, T_guess)[0]
        return T_bubble

    def dew_temperature(self, P, y, T_guess=None):
        """
        Calcule la temperature de rosee

        Returns :
        --------
        T_dew : float (K)
        """
        if T_guess is None:
            T_guess = np.sum([y[i] * self.compounds[i].Tb for i in range(self.n_comp)])

        def equation(T):
            T_val = T[0] if isinstance(T, np.ndarray) else T
            K = self.K_values(T_val, P)
            # Avoid division by zero
            return np.sum(y / (K + 1e-10)) - 1.0

        T_dew = fsolve(equation, T_guess)[0]
        return T_dew

    def mixture_enthalpy_liquid(self, T, x):
        """ Enthalpie du melange liquide (J/mol)"""
        h = np.sum([x[i] * self.compounds[i].enthalpy_liquid(T) for i in range(self.n_comp)])
        return h

    def mixture_enthalpy_vapor(self, T, y):
        """ Enthalpie du melange vapeur (J/mol )"""
        h = np.sum([y[i] * self.compounds[i].enthalpy_vapor(T) for i in range(self.n_comp)])
        return h
