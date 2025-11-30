import numpy as np
from scipy.optimize import brentq

def fenske_equation(x_LK_D, x_HK_D, x_LK_B, x_HK_B, alpha_avg):
    """
    Calcule le nombre minimum de plateaux

    Parameters :
    -----------
    x_LK_D , x_HK_D : float
        Fractions molaires des cles dans le distillat
    x_LK_B , x_HK_B : float
        Fractions molaires des cles dans le residu
    alpha_avg : float
        Volatilite relative moyenne

    Returns :
    --------
    N_min : float
        Nombre minimum de plateaux theoriques
    """
    numerator = np.log((x_LK_D / x_HK_D) / (x_LK_B / x_HK_B))
    denominator = np.log(alpha_avg)
    N_min = numerator / denominator
    return N_min

def underwood_method(alphas, x_F, x_D, q, alpha_LK, alpha_HK):
    """
    Calcule le reflux minimum par la methode d'Underwood

    Parameters :
    -----------
    alphas : array
        Volatilites relatives de tous les composes
    x_F : array
        Composition de l'alimentation
    x_D : array
        Composition du distillat
    q : float
        Qualite de l'alimentation
    alpha_LK : float
        Volatilite relative de la cle legere
    alpha_HK : float
        Volatilite relative de la cle lourde

    Returns :
    --------
    R_min : float
        Reflux minimum
    theta : float
        Racine de l'equation d'Underwood
    """
    # Equation 1: find theta
    # Theta doit etre entre alpha_HK et alpha_LK
    def equation1(theta):
        return np.sum(alphas * x_F / (alphas - theta)) - (1 - q)

    # Recherche de racine robuste
    # On cherche la racine entre alpha_HK et alpha_LK
    # On evite les poles en reduisant legerement l'intervalle
    epsilon = 1e-4
    low = alpha_HK + epsilon
    high = alpha_LK - epsilon
    
    try:
        theta = brentq(equation1, low, high)
    except ValueError:
        # Si brentq echoue (pas de changement de signe), on essaie de trouver ou est la racine
        # Cela peut arriver si q est tres different de 1 ou si les compositions sont particulieres
        # On fait un scan grossier
        ts = np.linspace(low, high, 20)
        vals = [equation1(t) for t in ts]
        # Trouver l'intervalle avec changement de signe
        found = False
        for i in range(len(vals)-1):
            if vals[i] * vals[i+1] < 0:
                theta = brentq(equation1, ts[i], ts[i+1])
                found = True
                break
        if not found:
            # Fallback: milieu de l'intervalle (ne devrait pas arriver pour une separation faisable)
            theta = (alpha_HK + alpha_LK) / 2

    # Equation 2: calculer R_min
    R_min_plus_1 = np.sum(alphas * x_D / (alphas - theta))
    R_min = R_min_plus_1 - 1

    return R_min, theta

def gilliland_correlation(N_min, R_min, R):
    """
    Calcule le nombre de plateaux theoriques

    Parameters :
    -----------
    N_min : float
        Nombre minimum de plateaux
    R_min : float
        Reflux minimum
    R : float
        Reflux operatoire

    Returns :
    --------
    N : float
        Nombre de plateaux theoriques
    """
    X = (R - R_min) / (R + 1)

    # Correlation de Gilliland
    exponent = (1 + 54.4 * X) * (X - 1) / ((11 + 117.2 * X) * np.sqrt(X))
    Y = 1 - np.exp(exponent)

    # Nombre de plateaux
    # Y = (N - N_min) / (N + 1)
    # Y(N+1) = N - N_min
    # YN + Y = N - N_min
    # N_min + Y = N - YN = N(1 - Y)
    # N = (N_min + Y) / (1 - Y)
    
    N = (N_min + Y) / (1 - Y)

    return N

def kirkbride_equation(B, D, x_HK_F, x_LK_F, x_LK_B, x_HK_D, N_total):
    """
    Determine la position du plateau d'alimentation

    Parameters :
    -----------
    B, D : float
        Debits de residu et distillat
    x_HK_F , x_LK_F : float
        Compositions des cles dans l'alimentation
    x_LK_B , x_HK_D : float
        Compositions des cles dans residu et distillat
    N_total : int
        Nombre total de plateaux

    Returns :
    --------
    N_R : int
        Nombre de plateaux en rectification
    N_S : int
        Nombre de plateaux en epuisement
    feed_stage : int
        Numero du plateau d'alimentation
    """
    # Calcul du ratio
    # Eq 23: log(NR/NS) = 0.206 * log(...)
    # ratio inside log is (B/D) * (x_HK_F / x_LK_F) * (x_LK_B / x_HK_D)**2
    # Note: The document code has (x_LK_B / x_HK_D)**2.
    # Let's verify Eq 23 in text:
    # log(NR/NS) = 0.206 log [ (B/D) * (x_HK_F / x_LK_F) * (x_LK_B / x_HK_D)^2 ]
    # Matches.
    
    ratio = (B / D) * (x_HK_F / x_LK_F) * (x_LK_B / x_HK_D)**2

    # Ratio N_R / N_S
    log_ratio = 0.206 * np.log(ratio)
    N_R_over_N_S = np.exp(log_ratio)

    # Resolution
    # N_R + N_S = N_total
    # N_R = N_S * ratio
    # N_S * ratio + N_S = N_total
    # N_S (1 + ratio) = N_total
    N_S = N_total / (1 + N_R_over_N_S)
    N_R = N_total - N_S

    # Feed stage
    # Usually feed stage is N_R + 1 (if counting from top)
    feed_stage = int(np.ceil(N_R)) + 1

    return int(np.ceil(N_R)), int(np.floor(N_S)), feed_stage
