import numpy as np
import plotly.graph_objects as go
from methode_matricielle.matrix_solver import MatrixDistillationSolver
import streamlit as st

def exemple_ethanol_eau():
    """
    Exemple de calcul pour un système éthanol-eau
    """
    # Création du solveur
    solver = MatrixDistillationSolver()

    # Configuration du problème
    components = ["Ethanol", "Eau"]
    n_stages = 12  # Nombre d'étages typique pour la séparation éthanol-eau
    feed_stage = 6  # Alimentation au milieu de la colonne
    
    # Volatilités relatives (à 78°C, pression atmosphérique)
    alpha_values = {
        "Ethanol": 2.25,  # Plus volatile
        "Eau": 1.0       # Composé de référence
    }
    
    # Composition de l'alimentation (mélange typique à 10% mol éthanol)
    feed_composition = {
        "Ethanol": 0.10,
        "Eau": 0.90
    }
    
    # Configuration du problème
    solver.setup_problem(components, n_stages, feed_stage, alpha_values)
    
    # Paramètres opératoires
    reflux_ratio = 3.0    # Taux de reflux
    boilup_ratio = 2.5    # Taux de rebouillage
    
    # Construction et analyse des matrices
    A, b = solver.build_matrix(reflux_ratio, boilup_ratio, feed_composition)
    
    # Calcul des profils
    profiles = solver.solve_profiles(reflux_ratio, boilup_ratio, feed_composition)
    
    # Affichage des résultats
    st.title("Exemple de Distillation Éthanol-Eau")
    
    # Informations sur le système
    st.header("Configuration du Système")
    st.write(f"Nombre d'étages: {n_stages}")
    st.write(f"Étage d'alimentation: {feed_stage}")
    st.write(f"Taux de reflux: {reflux_ratio}")
    st.write(f"Taux de rebouillage: {boilup_ratio}")
    
    # Composition d'alimentation
    st.subheader("Composition d'alimentation")
    for comp, value in feed_composition.items():
        st.write(f"{comp}: {value:.3f}")
    
    # Analyse de la matrice
    st.header("Analyse de la Matrice")
    st.text(solver.get_matrix_details())
    
    # Visualisation des matrices
    st.header("Visualisation des Matrices")
    fig_matrix = solver.visualize_matrices()
    st.plotly_chart(fig_matrix)
    
    # Profils de concentration
    st.header("Profils de Concentration")
    fig_profiles = solver.plot_composition_profiles(profiles)
    st.plotly_chart(fig_profiles)
    
    # Résultats de la séparation
    st.header("Résultats de la Séparation")
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Composition du Distillat")
        for comp in components:
            st.write(f"{comp}: {profiles[comp][0]:.4f}")
            
    with col2:
        st.subheader("Composition du Résidu")
        for comp in components:
            st.write(f"{comp}: {profiles[comp][-1]:.4f}")
    
    # Calcul du reflux minimum
    min_reflux = solver.calculate_minimum_reflux(feed_composition, components)
    st.write(f"\nReflux minimum calculé: {min_reflux:.2f}")
    st.write(f"Ratio reflux/reflux minimum: {reflux_ratio/min_reflux:.2f}")

if __name__ == "__main__":
    exemple_ethanol_eau()
