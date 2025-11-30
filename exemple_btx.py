"""
Exemple complet: Separation BTX (Benzene-Toluene-Xylene)

Ce script demonstre l'utilisation des methodes simplifiees et de la methode MESH
pour dimensionner une colonne de distillation multicomposants.
"""

import numpy as np
from compound import Compound
from thermodynamics import ThermodynamicPackage
from shortcut_methods import fenske_equation, underwood_method, gilliland_correlation, kirkbride_equation
from mesh_solver import MESHSolver
from visualizer import DistillationVisualizer

def main():
    print("="*70)
    print("SEPARATION BTX - DISTILLATION MULTICOMPOSANTS")
    print("="*70)
    print()

    # ========== DEFINITION DU SYSTEME ==========
    print("1. DEFINITION DU SYSTEME")
    print("-" * 70)

    # Composes
    compound_names = ['Benzene', 'Toluene', 'o-Xylene']
    compounds = [Compound(name) for name in compound_names]
    thermo = ThermodynamicPackage(compounds)

    # Conditions operatoires
    P = 101325  # Pa (1 atm)
    F = 100 / 3600  # kmol/h -> mol/s
    z_F = np.array([0.333, 0.333, 0.334])  # Composition alimentation

    print(f"Alimentation: {F * 3600:.1f} kmol/h")
    print(f"Pression: {P / 1e5:.3f} bar")
    print(f"Composition:")
    for i, name in enumerate(compound_names):
        print(f"  - {name}: {z_F[i] * 100:.1f}%")
    print()

    # Proprietes des composes
    print("Proprietes des composes:")
    print(f"{'Compose':<15} {'Tb (°C)':<10} {'Tc (K)':<10} {'Pc (bar)':<10}")
    for comp in compounds:
        print(f"{comp.name:<15} {comp.Tb - 273.15:<10.1f} {comp.Tc:<10.1f} {comp.Pc / 1e5:<10.2f}")
    print()

    # Calcul des volatilites relatives
    T_avg = np.mean([comp.Tb for comp in compounds])
    K_values = thermo.K_values(T_avg, P)
    alphas = K_values / K_values[1]  # Reference: Toluene

    print("Volatilites relatives (ref: Toluene):")
    for i, name in enumerate(compound_names):
        print(f"  - {name}: {alphas[i]:.2f}")
    print()

    # ========== SPECIFICATIONS DE SEPARATION ==========
    print("2. SPECIFICATIONS DE SEPARATION")
    print("-" * 70)

    # Objectif: 95% benzene dans distillat, 95% toluene+xylene dans residu
    # Cles: Benzene (LK), Toluene (HK)
    recovery_LK_D = 0.95  # 95% benzene dans distillat
    recovery_HK_B = 0.95  # 95% toluene dans residu

    # Bilans matieres
    F_benzene = F * z_F[0]
    F_toluene = F * z_F[1]
    F_xylene = F * z_F[2]

    D_benzene = F_benzene * recovery_LK_D
    B_benzene = F_benzene - D_benzene

    B_toluene = F_toluene * recovery_HK_B
    D_toluene = F_toluene - B_toluene

    # Xylene va presque tout dans le residu
    B_xylene = F_xylene * 0.999
    D_xylene = F_xylene - B_xylene

    D = D_benzene + D_toluene + D_xylene
    B = F - D

    x_D = np.array([D_benzene, D_toluene, D_xylene]) / D
    x_B = np.array([B_benzene, B_toluene, B_xylene]) / B

    print(f"Debit distillat: {D * 3600:.2f} kmol/h")
    print(f"Debit residu: {B * 3600:.2f} kmol/h")
    print()
    print("Composition distillat:")
    for i, name in enumerate(compound_names):
        print(f"  - {name}: {x_D[i] * 100:.2f}%")
    print()
    print("Composition residu:")
    for i, name in enumerate(compound_names):
        print(f"  - {name}: {x_B[i] * 100:.2f}%")
    print()

    # ========== METHODES SIMPLIFIEES ==========
    print("3. METHODES SIMPLIFIEES (SHORT-CUT)")
    print("-" * 70)

    # Methode de Fenske (Nmin)
    alpha_LK_HK = alphas[0] / alphas[1]  # Benzene/Toluene
    N_min = fenske_equation(x_D[0], x_D[1], x_B[0], x_B[1], alpha_LK_HK)
    print(f"Methode de Fenske:")
    print(f"  - Volatilite relative moyenne (Benzene/Toluene): {alpha_LK_HK:.2f}")
    print(f"  - Nombre minimum de plateaux (Nmin): {N_min:.1f}")
    print()

    # Methode d'Underwood (Rmin)
    q = 1.0  # Alimentation liquide saturee
    R_min, theta = underwood_method(alphas, z_F, x_D, q)
    print(f"Methode d'Underwood:")
    print(f"  - Reflux minimum (Rmin): {R_min:.2f}")
    print(f"  - Theta: {theta:.2f}")
    print()

    # Correlation de Gilliland (N)
    R_op = 1.3 * R_min  # Reflux operatoire = 1.3 * Rmin
    N_theoretical = gilliland_correlation(N_min, R_min, R_op)
    efficiency = 0.70
    N_actual = int(np.ceil(N_theoretical / efficiency))

    print(f"Correlation de Gilliland:")
    print(f"  - Reflux operatoire (1.3 * Rmin): {R_op:.2f}")
    print(f"  - Nombre de plateaux theoriques: {N_theoretical:.1f}")
    print(f"  - Efficacite des plateaux: {efficiency * 100:.0f}%")
    print(f"  - Nombre de plateaux reels: {N_actual}")
    print()

    # Equation de Kirkbride (position alimentation)
    N_R, N_S, feed_stage = kirkbride_equation(
        B, D, x_B[1], x_D[0], x_B[0], x_D[1], N_actual
    )
    print(f"Equation de Kirkbride:")
    print(f"  - Plateaux en rectification: {N_R}")
    print(f"  - Plateaux en epuisement: {N_S}")
    print(f"  - Plateau d'alimentation: {feed_stage}")
    print()

    # ========== SIMULATION RIGOUREUSE (MESH) ==========
    print("4. SIMULATION RIGOUREUSE (MESH)")
    print("-" * 70)

    # Initialisation du solveur MESH
    mesh_solver = MESHSolver(thermo, N_actual, feed_stage, P)
    mesh_solver.initialize(F, z_F, D, R_op)

    # Resolution
    print("Resolution du systeme MESH...")
    converged = mesh_solver.solve(max_iter=50, tol=1e-4)
    print()

    if converged:
        results = mesh_solver.get_results()

        print("RESULTATS DE LA SIMULATION:")
        print("-" * 70)
        print(f"Temperature en tete: {results['temperatures'][0] - 273.15:.1f} °C")
        print(f"Temperature en fond: {results['temperatures'][-1] - 273.15:.1f} °C")
        print()

        print("Composition distillat (simulation):")
        for i, name in enumerate(compound_names):
            print(f"  - {name}: {results['distillate']['composition'][i] * 100:.2f}%")
        print()

        print("Composition residu (simulation):")
        for i, name in enumerate(compound_names):
            print(f"  - {name}: {results['bottoms']['composition'][i] * 100:.2f}%")
        print()

        # ========== VISUALISATION ==========
        print("5. VISUALISATION DES RESULTATS")
        print("-" * 70)

        visualizer = DistillationVisualizer()

        # Profils de composition
        fig_comp = visualizer.plot_composition_profiles(
            results['stages'],
            results['compositions'],
            compound_names
        )
        fig_comp.write_html("btx_composition_profiles.html")
        print("Profils de composition sauvegardes: btx_composition_profiles.html")

        # Profil de temperature
        fig_temp = visualizer.plot_temperature_profile(
            results['stages'],
            results['temperatures']
        )
        fig_temp.write_html("btx_temperature_profile.html")
        print("Profil de temperature sauvegarde: btx_temperature_profile.html")
        print()

        # ========== TABLEAU RECAPITULATIF ==========
        print("6. TABLEAU RECAPITULATIF")
        print("-" * 70)
        print(f"{'Parametre':<30} {'Valeur':<15} {'Unite':<10}")
        print("-" * 70)
        print(f"{'Volatilite Benzene/Toluene':<30} {alpha_LK_HK:<15.2f} {'-':<10}")
        print(f"{'Nmin (Fenske)':<30} {N_min:<15.1f} {'-':<10}")
        print(f"{'Rmin (Underwood)':<30} {R_min:<15.2f} {'-':<10}")
        print(f"{'Rop (1.3 x Rmin)':<30} {R_op:<15.2f} {'-':<10}")
        print(f"{'Ntheorique (Gilliland)':<30} {N_theoretical:<15.1f} {'-':<10}")
        print(f"{'Nreel (E=0.70)':<30} {N_actual:<15} {'plateaux':<10}")
        print(f"{'Plateau alimentation':<30} {feed_stage:<15} {'-':<10}")
        print(f"{'Debit distillat':<30} {D * 3600:<15.2f} {'kmol/h':<10}")
        print(f"{'Debit residu':<30} {B * 3600:<15.2f} {'kmol/h':<10}")
        print(f"{'T tete':<30} {results['temperatures'][0] - 273.15:<15.1f} {'°C':<10}")
        print(f"{'T fond':<30} {results['temperatures'][-1] - 273.15:<15.1f} {'°C':<10}")
        print()

    print("="*70)
    print("SIMULATION TERMINEE")
    print("="*70)

if __name__ == "__main__":
    main()
