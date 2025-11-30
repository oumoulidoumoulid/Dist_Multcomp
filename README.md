# Distillation Multicomposants - Simulation Python

Ce projet implémente la modélisation et simulation de colonnes de distillation multicomposants en Python, conformément au cours "Modélisation et Simulation des Procédés" (PIC - UH1).

## 📋 Description

Le projet comprend:
- **Méthodes simplifiées** (Short-Cut Methods): Fenske, Underwood, Gilliland, Kirkbride
- **Méthode rigoureuse**: Système d'équations MESH (Material, Equilibrium, Summation, Heat)
- **Package thermodynamique**: Calculs d'équilibre liquide-vapeur avec la bibliothèque `thermo`
- **Visualisations interactives**: Profils de composition et température avec Plotly
- **Application web**: Interface Streamlit pour simulations interactives

## 🚀 Installation

### Prérequis
- Python 3.8 ou supérieur
- pip

### Installation des dépendances

```bash
pip install -r requirements.txt
```

## 📁 Structure du Projet

```
dist/
├── compound.py              # Classe Compound (propriétés des composés)
├── thermodynamics.py        # Package thermodynamique
├── shortcut_methods.py      # Méthodes simplifiées (Fenske, Underwood, etc.)
├── mesh_solver.py           # Solveur MESH rigoureux
├── visualizer.py            # Visualisations Plotly
├── exemple_btx.py           # Exemple complet: séparation BTX
├── app.py                   # Application Streamlit
├── requirements.txt         # Dépendances Python
└── README.md               # Ce fichier
```

## 🎯 Utilisation

### 1. Exemple BTX (Benzène-Toluène-Xylène)

Exécuter la simulation complète du système BTX:

```bash
python exemple_btx.py
```

Ce script:
- Calcule les paramètres de dimensionnement (Nmin, Rmin, N, position alimentation)
- Simule la colonne avec la méthode MESH
- Génère des visualisations HTML interactives

### 2. Application Streamlit Interactive

Lancer l'application web:

```bash
streamlit run app.py
```

L'application permet de:
- Sélectionner les composés du mélange
- Définir les conditions opératoires (pression, débit, composition)
- Spécifier les clés de séparation et récupérations
- Calculer automatiquement les paramètres de design
- Exporter un rapport de dimensionnement

## 📊 Méthodes Implémentées

### Méthodes Simplifiées

1. **Fenske**: Calcul du nombre minimum de plateaux (Nmin)
   ```
   Nmin = log[(xLK/xHK)D / (xLK/xHK)B] / log(αLK,HK)
   ```

2. **Underwood**: Calcul du reflux minimum (Rmin)
   - Résolution de l'équation pour θ
   - Calcul de Rmin à partir de θ

3. **Gilliland**: Corrélation empirique pour le nombre de plateaux
   ```
   Y = f(X) où X = (R - Rmin)/(R + 1)
   ```

4. **Kirkbride**: Position optimale du plateau d'alimentation
   ```
   log(NR/NS) = 0.206 log[(B/D)(xHK,F/xLK,F)(xLK,B/xHK,D)²]
   ```

### Méthode Rigoureuse MESH

Résolution simultanée des équations:
- **M**: Bilans matières pour chaque composé sur chaque plateau
- **E**: Équilibre liquide-vapeur (yi = Ki·xi)
- **S**: Contraintes de sommation (Σxi = 1, Σyi = 1)
- **H**: Bilans énergétiques

## 📈 Visualisations

Le projet génère des graphiques interactifs Plotly:

1. **Profils de composition**: Fractions molaires liquide et vapeur vs numéro de plateau
2. **Profil de température**: Température vs numéro de plateau

Les graphiques sont sauvegardés en HTML et peuvent être ouverts dans un navigateur.

## 🧪 Exemple: Séparation BTX

**Système**:
- Alimentation: 100 kmol/h (33.3% Benzène, 33.3% Toluène, 33.4% Xylène)
- Pression: 1.013 bar
- Objectif: 95% Benzène dans distillat, 95% Toluène dans résidu

**Résultats typiques**:
- Nmin: ~7 plateaux
- Rmin: ~1.4
- Nréel: ~22 plateaux (efficacité 70%)
- Plateau alimentation: ~11
- T tête: ~80°C
- T fond: ~139°C

## 📚 Références

- Cours: "Modélisation et Simulation des Procédés" - Prof. BAKHER Zine Elabidine
- Filière: Procédés et Ingénierie Chimique (PIC)
- Université: UH1
- Année: 2024-2025

## 🔧 Dépendances

- `thermo`: Propriétés thermodynamiques et équilibre de phases
- `scipy`: Résolution d'équations non-linéaires
- `numpy`: Calculs numériques
- `plotly`: Visualisations interactives
- `pandas`: Manipulation de données
- `streamlit`: Application web interactive

## 📝 Notes

- Les méthodes simplifiées supposent des volatilités relatives constantes et des mélanges idéaux
- La méthode MESH est plus précise mais nécessite plus de temps de calcul
- L'application Streamlit permet une exploration interactive des paramètres

## 🎓 Compétences Développées

- Modélisation thermodynamique de systèmes multicomposants
- Algorithmes de résolution numérique (Newton-Raphson, méthodes tridiagonales)
- Programmation orientée objet en Python
- Visualisation scientifique interactive
- Développement d'applications web pour l'ingénierie

## 📧 Contact

Pour toute question concernant ce projet, veuillez contacter le département de Procédés et Ingénierie Chimique (PIC) - UH1.
