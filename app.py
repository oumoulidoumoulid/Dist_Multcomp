"""
Application Streamlit pour Distillation Multicomposants
ALIGNEE AVEC LES RESULTATS DU PDF DE REFERENCE

Resultats attendus pour BTX (Table 2 du PDF):
- α benzene/toluene: 2.41
- Nmin: 6.8
- Rmin: 1.85
- Rop: 2.41 (1.3 × Rmin)
- N theorique: 13.2
- N reel: 19 plateaux (E=70%)
- Plateau alimentation: 10
- T tete: 80°C
- T fond: 140°C
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from compound import Compound
from thermodynamics import ThermodynamicPackage
from shortcut_methods import fenske_equation, underwood_method, gilliland_correlation, kirkbride_equation

st.set_page_config(page_title="Distillation Multicomposants", layout="wide", page_icon="🧪")

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.5rem;
        font-weight: bold;
        color: #2ca02c;
        margin-top: 2rem;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">🧪 Simulation de Distillation Multicomposants</div>', unsafe_allow_html=True)

# ========== SIDEBAR ==========
st.sidebar.header("⚙️ Parametres de Simulation")

# Composes
st.sidebar.subheader("Composes")
available_compounds = ['Benzene', 'Toluene', 'o-Xylene', 'Ethanol', 'Water', 'Methanol', 'Acetone']

compound_names = st.sidebar.multiselect(
    "Selectionner les composes",
    available_compounds,
    default=['Benzene', 'Toluene', 'o-Xylene']
)

if len(compound_names) < 2:
    st.warning("⚠️ Veuillez selectionner au moins 2 composes")
    st.stop()

n_comp = len(compound_names)

# Conditions operatoires
st.sidebar.subheader("Conditions Operatoires")
P_top_bar = st.sidebar.number_input("Pression Tete (bar)", value=1.013, min_value=0.1, max_value=10.0, step=0.001, format="%.3f")
P_drop_bar = st.sidebar.number_input("Perte de charge (bar)", value=0.0, min_value=0.0, max_value=2.0, step=0.01)
P_top = P_top_bar * 1e5

F_kmolh = st.sidebar.number_input("Debit alimentation (kmol/h)", value=100.0, min_value=1.0, step=1.0)
F = F_kmolh / 3600

q = st.sidebar.slider("Qualite Alimentation (q)", 0.0, 1.0, 1.0, 0.01, help="1=Liquide sature, 0=Vapeur saturee")

z_F = []
for i, name in enumerate(compound_names):
    default_val = 100.0 / n_comp
    z_i_pct = st.sidebar.number_input(
        f"{name}",
        value=default_val,
        min_value=0.0,
        max_value=100.0,
        step=0.1,
        key=f"z_{i}",
        format="%.1f"
    )
    z_F.append(z_i_pct / 100.0)

z_F = np.array(z_F)
sum_z = np.sum(z_F)
if abs(sum_z - 1.0) > 0.01:
    st.sidebar.warning(f"⚠️ Somme = {sum_z*100:.1f}% (normalisation auto)")
    z_F = z_F / sum_z

# Specifications
st.sidebar.subheader("Specifications")
LK_idx = st.sidebar.selectbox("Cle leger (LK)", range(n_comp), format_func=lambda x: compound_names[x])
HK_idx = st.sidebar.selectbox("Cle lourd (HK)", range(n_comp), format_func=lambda x: compound_names[x], index=min(1, n_comp-1))

recovery_LK = st.sidebar.slider("Recuperation LK dans distillat (%)", 80, 99, 95) / 100
recovery_HK = st.sidebar.slider("Recuperation HK dans residu (%)", 80, 99, 95) / 100

# Design
st.sidebar.subheader("Design")
R_factor = st.sidebar.slider("Facteur de reflux (x Rmin)", 1.1, 3.0, 1.3, 0.1)
efficiency = st.sidebar.slider("Efficacite plateaux (%)", 50, 100, 70) / 100

# ========== SIMULATION ==========
if st.sidebar.button("🚀 Lancer la Simulation", type="primary"):
    
    with st.spinner("Calcul en cours..."):
        try:
            compounds = [Compound(name) for name in compound_names]
            thermo = ThermodynamicPackage(compounds)
        except Exception as e:
            st.error(f"❌ Erreur: {e}")
            st.stop()

    # ========== PROPRIETES ==========
    st.markdown('<div class="section-header">📊 Proprietes des Composes</div>', unsafe_allow_html=True)
    
    props_data = []
    for comp in compounds:
        props_data.append({
            'Compose': comp.name,
            'Tb (°C)': f"{comp.Tb - 273.15:.1f}",
            'Tc (K)': f"{comp.Tc:.1f}",
            'Pc (bar)': f"{comp.Pc / 1e5:.2f}",
            'MW (g/mol)': f"{comp.MW:.2f}"
        })
    
    df_props = pd.DataFrame(props_data)
    st.dataframe(df_props, use_container_width=True, hide_index=True)

    # Volatilites (reference: composant le plus lourd)
    T_ref = 383.15  # 110°C - temperature de reference pour BTX
    K_values = thermo.K_values(T_ref, P_top)
    
    # Trouver le composant de reference (le plus lourd = K le plus petit)
    ref_idx = np.argmin(K_values)
    alphas = K_values / K_values[ref_idx]

    st.subheader(f"Volatilites Relatives (ref: {compound_names[ref_idx]})")
    alpha_data = pd.DataFrame({
        'Compose': compound_names,
        'α': [f"{a:.2f}" for a in alphas]
    })
    st.dataframe(alpha_data, use_container_width=True, hide_index=True)

    # ========== BILANS MATIERES ==========
    st.markdown('<div class="section-header">⚖️ Bilans Matieres</div>', unsafe_allow_html=True)

    F_comp = F * z_F
    D_comp = np.zeros(n_comp)
    B_comp = np.zeros(n_comp)

    # Cles
    D_comp[LK_idx] = F_comp[LK_idx] * recovery_LK
    B_comp[LK_idx] = F_comp[LK_idx] - D_comp[LK_idx]
    B_comp[HK_idx] = F_comp[HK_idx] * recovery_HK
    D_comp[HK_idx] = F_comp[HK_idx] - B_comp[HK_idx]

    # Non-cles
    for i in range(n_comp):
        if i != LK_idx and i != HK_idx:
            if alphas[i] > alphas[LK_idx]:
                D_comp[i] = F_comp[i] * 0.999
                B_comp[i] = F_comp[i] - D_comp[i]
            else:
                B_comp[i] = F_comp[i] * 0.999
                D_comp[i] = F_comp[i] - B_comp[i]

    D = np.sum(D_comp)
    B = np.sum(B_comp)
    x_D = D_comp / D
    x_B = B_comp / B

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Debit Alimentation", f"{F * 3600:.1f} kmol/h")
    with col2:
        st.metric("Debit Distillat", f"{D * 3600:.1f} kmol/h")
    with col3:
        st.metric("Debit Residu", f"{B * 3600:.1f} kmol/h")

    balance_data = []
    for i, name in enumerate(compound_names):
        recup_D = (D_comp[i] / F_comp[i] * 100) if F_comp[i] > 0 else 0
        recup_B = (B_comp[i] / F_comp[i] * 100) if F_comp[i] > 0 else 0
        balance_data.append({
            'Compose': name,
            'Alim. (kmol/h)': f"{F_comp[i] * 3600:.1f}",
            'Dist. (kmol/h)': f"{D_comp[i] * 3600:.1f}",
            'Residu (kmol/h)': f"{B_comp[i] * 3600:.1f}",
            'Recup. D (%)': f"{recup_D:.1f}",
            'Recup. B (%)': f"{recup_B:.1f}"
        })
    
    df_balance = pd.DataFrame(balance_data)
    st.dataframe(df_balance, use_container_width=True, hide_index=True)

    # ========== METHODES SIMPLIFIEES ==========
    st.markdown('<div class="section-header">📐 Dimensionnement (Methodes Simplifiees)</div>', unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🔹 Fenske")
        alpha_LK_HK = alphas[LK_idx] / alphas[HK_idx]
        N_min = fenske_equation(x_D[LK_idx], x_D[HK_idx], x_B[LK_idx], x_B[HK_idx], alpha_LK_HK)
        st.metric("α (LK/HK)", f"{alpha_LK_HK:.2f}")
        st.metric("Nmin", f"{N_min:.1f} plateaux")

    with col2:
        st.subheader("🔹 Underwood")
        try:
            # Passer les volatilites des cles explicitement
            alpha_LK_val = alphas[LK_idx]
            alpha_HK_val = alphas[HK_idx]
            R_min, theta = underwood_method(alphas, z_F, x_D, q, alpha_LK_val, alpha_HK_val)
            st.metric("Rmin", f"{R_min:.2f}")
            st.metric("θ", f"{theta:.2f}")
        except Exception as e:
            st.error(f"Erreur: {e}")
            R_min = 1.0
            theta = 1.0

    col1, col2, col3 = st.columns(3)
    
    R_op = R_factor * R_min
    N_theoretical = gilliland_correlation(N_min, R_min, R_op)
    N_actual = int(np.ceil(N_theoretical / efficiency))
    
    with col1:
        st.metric("Reflux operatoire", f"{R_op:.2f}")
    with col2:
        st.metric("N theorique", f"{N_theoretical:.1f}")
    with col3:
        st.metric("N reel", f"{N_actual} plateaux")

    N_R, N_S, feed_stage = kirkbride_equation(B, D, x_B[HK_idx], x_D[LK_idx], x_B[LK_idx], x_D[HK_idx], N_actual)

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Plateaux rectification", f"{N_R}")
    with col2:
        st.metric("Plateaux epuisement", f"{N_S}")
    with col3:
        st.metric("Plateau alimentation", f"{feed_stage}")

    # ========== PROFILS ==========
    st.markdown('<div class="section-header">📈 Profils de Composition et Temperature</div>', unsafe_allow_html=True)
    
    st.info("💡 Profils calculés avec équilibres liquide-vapeur stage-by-stage (Heuristique améliorée)")
    
    stages = np.arange(1, N_actual + 1)
    x_profiles = np.zeros((N_actual, n_comp))
    y_profiles = np.zeros((N_actual, n_comp))
    T_profile = np.zeros(N_actual)
    P_profile = np.linspace(P_top, P_top + P_drop_bar * 1e5, N_actual)
    
    # Condenseur (Plateau 1)
    x_profiles[0, :] = x_D
    y_profiles[0, :] = x_D # Hypothese condenseur total pour affichage
    T_profile[0] = thermo.bubble_temperature(P_profile[0], x_D)
    
    # Rebouilleur (Plateau N)
    x_profiles[-1, :] = x_B
    T_profile[-1] = thermo.bubble_temperature(P_profile[-1], x_B)
    K_bottom = thermo.K_values(T_profile[-1], P_profile[-1])
    y_profiles[-1, :] = K_bottom * x_B
    y_profiles[-1, :] /= np.sum(y_profiles[-1, :])
    
    # Estimation composition plateau alimentation (z_F approx)
    x_feed = z_F
    
    # Plateaux intermediaires
    # Rectification: de 1 à feed_stage-1
    # Epuisement: de feed_stage à N-2
    
    for j in range(1, N_actual - 1):
        if j < feed_stage - 1:
            # Section Rectification (D -> F)
            # Fraction de progression (0 au condenseur, 1 au plateau alim)
            frac = j / (feed_stage - 1)
            
            for i in range(n_comp):
                # Interpolation non-lineaire entre x_D et x_feed
                # On utilise alpha pour courber le profil
                val_interp = x_D[i] * (1 - frac) + x_feed[i] * frac
                
                # Effet de separation
                if alphas[i] > 1.0:
                    # Composant leger: profil concave ou convexe selon position
                    x_profiles[j, i] = val_interp * (alphas[i] ** (0.2 * (1-frac)))
                else:
                    x_profiles[j, i] = val_interp * (alphas[i] ** (0.2 * frac))
                    
        else:
            # Section Epuisement (F -> B)
            # Fraction de progression (0 au plateau alim, 1 au rebouilleur)
            frac = (j - (feed_stage - 1)) / (N_actual - feed_stage)
            
            for i in range(n_comp):
                # Interpolation entre x_feed et x_B
                val_interp = x_feed[i] * (1 - frac) + x_B[i] * frac
                
                if alphas[i] > 1.0:
                    x_profiles[j, i] = val_interp * (alphas[i] ** (-0.2 * frac))
                else:
                    x_profiles[j, i] = val_interp * (alphas[i] ** (0.2 * (1-frac)))

        # Normalisation
        x_profiles[j, :] /= np.sum(x_profiles[j, :])
        
        # Equilibre
        T_profile[j] = thermo.bubble_temperature(P_profile[j], x_profiles[j, :])
        K_j = thermo.K_values(T_profile[j], P_profile[j])
        y_profiles[j, :] = K_j * x_profiles[j, :]
        y_profiles[j, :] /= np.sum(y_profiles[j, :])

    # Lissage au plateau d'alimentation pour eviter discontinuité
    # Moyenne mobile sur 3 points autour de l'alimentation
    if feed_stage > 1 and feed_stage < N_actual - 1:
        idx_feed = feed_stage - 1 # 0-indexed
        x_profiles[idx_feed, :] = (x_profiles[idx_feed-1, :] + x_profiles[idx_feed+1, :]) / 2
        x_profiles[idx_feed, :] /= np.sum(x_profiles[idx_feed, :])
        T_profile[idx_feed] = thermo.bubble_temperature(P_profile[idx_feed], x_profiles[idx_feed, :])
        K_feed = thermo.K_values(T_profile[idx_feed], P_profile[idx_feed])
        y_profiles[idx_feed, :] = K_feed * x_profiles[idx_feed, :]
        y_profiles[idx_feed, :] /= np.sum(y_profiles[idx_feed, :])

    # ========== GRAPHIQUES (STYLE MATPLOTLIB) ==========
    st.markdown('<div class="section-header">📈 Visualisations (Style Reference)</div>', unsafe_allow_html=True)
    
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    
    # Configuration du style
    plt.style.use('default')
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.3
    
    # 1. BILANS MATIERES
    st.subheader("1. Bilans Matières")
    fig1 = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.5])
    
    # Debits
    ax1 = plt.subplot(gs[0])
    flux = ['Alimentation', 'Distillat', 'Résidu']
    debits = [F*3600, D*3600, B*3600]
    colors_flux = ['#4444ff', '#44aa44', '#ff4444']
    
    bars = ax1.bar(flux, debits, color=colors_flux, edgecolor='black', alpha=0.8)
    ax1.set_ylabel('Débit (kmol/h)', fontweight='bold')
    ax1.set_title('Débits des flux', fontweight='bold')
    ax1.set_ylim(0, max(debits)*1.2)
    
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}\nkmol/h',
                ha='center', va='bottom', fontweight='bold')
    
    # Compositions
    ax2 = plt.subplot(gs[1])
    x_pos = np.arange(n_comp)
    width = 0.25
    
    rects1 = ax2.bar(x_pos - width, z_F, width, label='Alimentation', color='#4444ff', edgecolor='black', alpha=0.8)
    rects2 = ax2.bar(x_pos, x_D, width, label='Distillat', color='#44aa44', edgecolor='black', alpha=0.8)
    rects3 = ax2.bar(x_pos + width, x_B, width, label='Résidu', color='#ff4444', edgecolor='black', alpha=0.8)
    
    ax2.set_ylabel('Fraction molaire', fontweight='bold')
    ax2.set_title('Compositions des flux', fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(compound_names, rotation=15)
    ax2.legend()
    ax2.set_ylim(0, 1.0)
    
    st.pyplot(fig1)
    
    # 2. PROFILS DE COMPOSITION
    st.subheader("2. Profils de Composition")
    fig2, (ax_liq, ax_vap) = plt.subplots(1, 2, figsize=(14, 8))
    
    colors_comp = ['#8dd3c7', '#b3de69', '#fdb462', '#fb8072', '#80b1d3']
    
    # Liquide
    for i, name in enumerate(compound_names):
        ax_liq.plot(x_profiles[:, i], stages, 'o-', label=name, 
                   color=colors_comp[i % len(colors_comp)], linewidth=2, markersize=5)
    
    ax_liq.axhline(y=feed_stage, color='blue', linestyle='--', linewidth=2, label='Plateau alimentation')
    ax_liq.set_title('Phase Liquide', fontweight='bold')
    ax_liq.set_xlabel('Fraction molaire liquide (x)', fontweight='bold')
    ax_liq.set_ylabel('Numéro de plateau', fontweight='bold')
    ax_liq.invert_yaxis()
    ax_liq.legend()
    ax_liq.grid(True, alpha=0.3)
    
    # Vapeur
    for i, name in enumerate(compound_names):
        ax_vap.plot(y_profiles[:, i], stages, 's-', label=name,
                   color=colors_comp[i % len(colors_comp)], linewidth=2, markersize=5)
                   
    ax_vap.axhline(y=feed_stage, color='blue', linestyle='--', linewidth=2, label='Plateau alimentation')
    ax_vap.set_title('Phase Vapeur', fontweight='bold')
    ax_vap.set_xlabel('Fraction molaire vapeur (y)', fontweight='bold')
    ax_vap.set_ylabel('Numéro de plateau', fontweight='bold')
    ax_vap.invert_yaxis()
    ax_vap.legend()
    ax_vap.grid(True, alpha=0.3)
    
    st.pyplot(fig2)
    
    # 3. PROFIL DE TEMPERATURE
    st.subheader("3. Profil de Température")
    fig3, ax_temp = plt.subplots(figsize=(10, 8))
    
    ax_temp.plot(T_profile - 273.15, stages, 'o-', color='#ff4400', linewidth=3, markersize=8, label='Température')
    ax_temp.axhline(y=feed_stage, color='blue', linestyle='--', linewidth=2, label=f'Plateau alimentation ({feed_stage})')
    
    # Annotations
    ax_temp.text(T_profile[0]-273.15 + 1, 1, f"{T_profile[0]-273.15:.1f}°C", 
                color='#8B0000', fontweight='bold', va='center')
    ax_temp.text(T_profile[-1]-273.15 + 1, N_actual, f"{T_profile[-1]-273.15:.1f}°C", 
                color='#8B0000', fontweight='bold', va='center')
    
    ax_temp.set_title('Profil de Température dans la Colonne', fontweight='bold', fontsize=14)
    ax_temp.set_xlabel('Température (°C)', fontweight='bold')
    ax_temp.set_ylabel('Numéro de plateau', fontweight='bold')
    ax_temp.invert_yaxis()
    ax_temp.legend()
    ax_temp.grid(True, alpha=0.3)
    
    st.pyplot(fig3)
    
    # 4. EFFET DU REFLUX
    st.subheader("4. Effet du Rapport de Reflux")
    
    # Calcul de la courbe N vs R
    R_ratios = np.linspace(1.05, 3.0, 50)
    N_values = []
    
    for r_ratio in R_ratios:
        R_val = r_ratio * R_min
        N_val = gilliland_correlation(N_min, R_min, R_val)
        N_values.append(N_val)
        
    fig4, ax_reflux = plt.subplots(figsize=(12, 7))
    
    ax_reflux.plot(R_ratios, N_values, 'b-', linewidth=3, label='Courbe N vs R/R_min')
    
    # Point de fonctionnement
    ax_reflux.plot(R_factor, N_theoretical, 'ro', markersize=12, label=f'Point de fonctionnement\n(R={R_op:.2f}, N={N_theoretical:.1f})')
    
    # Asymptotes
    ax_reflux.axvline(x=R_factor, color='green', linestyle='--', linewidth=2, label=f'R = {R_factor}×R_min')
    ax_reflux.axhline(y=N_min, color='red', linestyle='--', linewidth=2, label=f'N_min = {N_min:.1f}')
    
    ax_reflux.set_title('Effet du rapport de reflux sur le nombre de plateaux', fontweight='bold', fontsize=14)
    ax_reflux.set_xlabel('R / R_min', fontweight='bold')
    ax_reflux.set_ylabel('Nombre de plateaux théoriques', fontweight='bold')
    ax_reflux.legend()
    ax_reflux.grid(True, alpha=0.3)
    
    st.pyplot(fig4)

    # ========== BILANS ENERGETIQUES ==========
    st.markdown('<div class="section-header">⚡ Bilans Energetiques</div>', unsafe_allow_html=True)
    
    T_condenser = T_profile[0]
    T_reboiler = T_profile[-1]
    
    lambda_condenser = 0.0
    for i, comp in enumerate(compounds):
        comp.chem.T = T_condenser
        lambda_i = comp.chem.Hvap if comp.chem.Hvap else 30000
        lambda_condenser += x_D[i] * lambda_i
    
    lambda_reboiler = 0.0
    for i, comp in enumerate(compounds):
        comp.chem.T = T_reboiler
        lambda_i = comp.chem.Hvap if comp.chem.Hvap else 30000
        lambda_reboiler += x_B[i] * lambda_i
    
    V_top = (R_op + 1) * D
    Q_condenser_kW = -(V_top * lambda_condenser) / 1000
    Q_reboiler_kW = (V_top * lambda_reboiler) / 1000
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("❄️ Condenseur")
        st.metric("Duty", f"{Q_condenser_kW:.2f} kW")
        st.metric("Temperature", f"{T_condenser - 273.15:.1f} °C")
        st.metric("Lambda moyen", f"{lambda_condenser/1000:.1f} kJ/mol")
    
    with col2:
        st.subheader("🔥 Rebouilleur")
        st.metric("Duty", f"{Q_reboiler_kW:.2f} kW")
        st.metric("Temperature", f"{T_reboiler - 273.15:.1f} °C")
        st.metric("Lambda moyen", f"{lambda_reboiler/1000:.1f} kJ/mol")

    # ========== TABLEAUX ==========
    st.markdown('<div class="section-header">📋 Tableau Recapitulatif</div>', unsafe_allow_html=True)
    
    summary_data = {
        'Parametre': [
            f'Volatilite α {compound_names[LK_idx]}/{compound_names[HK_idx]}',
            'Nmin (Fenske)',
            'Rmin (Underwood)',
            f'Rop ({R_factor}× Rmin)',
            'Ntheorique (Gilliland)',
            f'Nreel (E={efficiency*100:.0f}%)',
            'Plateau alimentation',
            'Debit distillat',
            'Debit residu',
            'Temperature tete',
            'Temperature fond',
            'Duty condenseur',
            'Duty rebouilleur'
        ],
        'Valeur': [
            f"{alpha_LK_HK:.2f}",
            f"{N_min:.1f}",
            f"{R_min:.2f}",
            f"{R_op:.2f}",
            f"{N_theoretical:.1f}",
            f"{N_actual}",
            f"{feed_stage}",
            f"{D * 3600:.1f} kmol/h",
            f"{B * 3600:.1f} kmol/h",
            f"{T_profile[0] - 273.15:.1f} °C",
            f"{T_profile[-1] - 273.15:.1f} °C",
            f"{Q_condenser_kW:.2f} kW",
            f"{Q_reboiler_kW:.2f} kW"
        ]
    }
    
    df_summary = pd.DataFrame(summary_data)
    st.dataframe(df_summary, use_container_width=True, hide_index=True)
    
    # Export
    csv_summary = df_summary.to_csv(index=False)
    st.download_button(
        label="📄 Telecharger Recapitulatif (CSV)",
        data=csv_summary,
        file_name="distillation_summary.csv",
        mime="text/csv"
    )

else:
    st.info("👈 Configurez les parametres et cliquez sur 'Lancer la Simulation'")
    
    st.markdown("### 📚 Exemple: Separation BTX (Valeurs de Reference)")
    st.markdown("""
    **Systeme**: Benzene-Toluene-o-Xylene
    - Alimentation: 100 kmol/h (33.3% chaque)
    - Pression: 1.013 bar
    - Objectif: 95% Benzene dans distillat
    
    **Resultats attendus (Table 2 du PDF)**:
    - α benzene/toluene: **2.41**
    - Nmin: **6.8**
    - Rmin: **1.85**
    - Rop: **2.41** (1.3 × Rmin)
    - N theorique: **13.2**
    - N reel: **19 plateaux** (E=70%)
    - Plateau alimentation: **10**
    - T tete: **80°C**
    - T fond: **140°C**
    """)
