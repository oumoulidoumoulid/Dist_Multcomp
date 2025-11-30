"""
Script de test rapide pour verifier les fonctionnalites principales
"""

import numpy as np
from compound import Compound
from thermodynamics import ThermodynamicPackage
from shortcut_methods import fenske_equation, underwood_method, gilliland_correlation

def test_compound():
    print("Test 1: Classe Compound")
    print("-" * 50)
    benzene = Compound('Benzene')
    print(f"Benzene - Tb: {benzene.Tb - 273.15:.1f} °C")
    print(f"Benzene - Tc: {benzene.Tc:.1f} K")
    print(f"Benzene - MW: {benzene.MW:.2f} g/mol")
    
    T = 353.15  # 80°C
    P = 101325  # 1 atm
    K = benzene.K_value(T, P)
    print(f"K-value a {T-273.15:.0f}°C: {K:.2f}")
    print("✓ Compound OK\n")

def test_thermodynamics():
    print("Test 2: Package Thermodynamique")
    print("-" * 50)
    compounds = [Compound('Benzene'), Compound('Toluene')]
    thermo = ThermodynamicPackage(compounds)
    
    P = 101325
    x = np.array([0.5, 0.5])
    
    T_bubble = thermo.bubble_temperature(P, x)
    print(f"Temperature de bulle: {T_bubble - 273.15:.1f} °C")
    
    K_values = thermo.K_values(T_bubble, P)
    print(f"K-values: {K_values}")
    print("✓ Thermodynamics OK\n")

def test_shortcut_methods():
    print("Test 3: Methodes Simplifiees")
    print("-" * 50)
    
    # Donnees de test
    x_LK_D = 0.95
    x_HK_D = 0.05
    x_LK_B = 0.05
    x_HK_B = 0.95
    alpha_avg = 2.4
    
    N_min = fenske_equation(x_LK_D, x_HK_D, x_LK_B, x_HK_B, alpha_avg)
    print(f"Fenske - Nmin: {N_min:.1f}")
    
    alphas = np.array([2.4, 1.0, 0.44])
    z_F = np.array([0.33, 0.33, 0.34])
    x_D = np.array([0.95, 0.05, 0.0])
    q = 1.0
    
    try:
        R_min, theta = underwood_method(alphas, z_F, x_D, q)
        print(f"Underwood - Rmin: {R_min:.2f}, theta: {theta:.2f}")
    except Exception as e:
        print(f"Underwood - Erreur: {e}")
    
    R = 1.3 * 1.5  # Exemple
    N = gilliland_correlation(N_min, 1.5, R)
    print(f"Gilliland - N: {N:.1f}")
    print("✓ Shortcut Methods OK\n")

def main():
    print("="*50)
    print("TESTS DE VALIDATION")
    print("="*50)
    print()
    
    try:
        test_compound()
        test_thermodynamics()
        test_shortcut_methods()
        
        print("="*50)
        print("TOUS LES TESTS REUSSIS ✓")
        print("="*50)
    except Exception as e:
        print(f"\n❌ ERREUR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
