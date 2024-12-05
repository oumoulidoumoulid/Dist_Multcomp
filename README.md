# Projet de Distillation Multicomposée

Ce projet implémente deux méthodes de calcul pour la distillation multicomposée :
1. Méthode FUG (Fugacité)
2. Méthode Matricielle

## Structure du Projet

```
distillation_project/
├── methode_fug/
│   ├── __init__.py
│   ├── fugacity_calculator.py
│   ├── thermodynamic_properties.py
│   └── utils.py
├── methode_matricielle/
│   ├── __init__.py
│   ├── matrix_solver.py
│   ├── equilibrium_calculator.py
│   └── utils.py
├── interface/
│   ├── app.py
│   └── styles/
├── tests/
└── data/

## Installation

1. Créer un environnement virtuel Python :
```bash
python -m venv venv
source venv/bin/activate  # Sur Windows: venv\Scripts\activate
```

2. Installer les dépendances :
```bash
pip install -r requirements.txt
```

## Utilisation

1. Pour la méthode FUG :
```bash
python -m interface.app --method fug
```

2. Pour la méthode Matricielle :
```bash
python -m interface.app --method matrix
```

## Fonctionnalités

### Méthode FUG
- Calcul des coefficients de fugacité
- Équilibres liquide-vapeur
- Visualisation des résultats
- Optimisation des paramètres

### Méthode Matricielle
- Résolution matricielle des équations de bilan
- Calcul des profils de concentration
- Visualisation des résultats
- Analyse de sensibilité
