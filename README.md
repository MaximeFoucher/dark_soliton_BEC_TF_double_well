# 🌊 Collision de Solitons Sombres dans un Condensat de Bose-Einstein Binaire

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![NumPy](https://img.shields.io/badge/NumPy-1.20%2B-orange)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.3%2B-green)](https://matplotlib.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

> **Simulation numérique de la collision entre deux solitons sombres d'espèces atomiques distinctes (Rb-87 et Cs-133) dans un piège harmonique quasi-1D.**

---

## 📋 Table des matières

- [Description](#-description)
- [Physique du système](#-physique-du-système)
- [Installation](#-installation)
- [Utilisation](#-utilisation)
- [Structure du code](#-structure-du-code)
- [Paramètres physiques](#-paramètres-physiques)
- [Résultats](#-résultats)
- [Références](#-références)
- [Auteurs](#-auteurs)

---

## 📖 Description

Ce projet simule la **dynamique de collision entre deux solitons sombres** portés par des atomes de rubidium, potassium, sodium et de césium dans un condensat de Bose-Einstein (BEC) binaire en régime quasi-1D. 

### Objectifs scientifiques

- Étudier l'influence du **paramètre de miscibilité** Δ = g₁₂² − g₁₁g₂₂ sur la collision
- Identifier les régimes de comportement (miscible, immiscible, intégrable)
- Valider numériquement les prédictions théoriques sur les mélanges quantiques

### Méthode numérique

Résolution des **équations de Gross-Pitaevskii couplées** via la méthode **split-step spectrale** (Beam Propagation Method).

---

## ⚛️ Physique du système

### Équations de Gross-Pitaevskii couplées

Le système est décrit par deux GPE couplées :
```
iℏ ∂ψ₁/∂t = [−ℏ²/2m₁ ∇² + V_ext + g₁₁|ψ₁|² + g₁₂|ψ₂|²] ψ₁
iℏ ∂ψ₂/∂t = [−ℏ²/2m₂ ∇² + V_ext + g₂₂|ψ₂|² + g₂₁|ψ₁|²] ψ₂
```

### Couplages inter-atomiques

Les constantes d'interaction sont calculées via :

1. **Couplage 3D** (formule de Dalibard) :
```
   g₃D = (4π ℏ² aₛ) / m
```

2. **Réduction quasi-1D** (formule d'Olshanii) :
```
   g₁D = 2 ℏ ω_⊥ aₛ
```

3. **Adimensionnement** :
```
   g_adim = g₁D / (ℏ ω_z a_z)
```

### Critère de miscibilité

**Δ = g₁₂² − g₁₁ × g₂₂**

- **Δ < 0** : Régime **miscible** (mélange homogène)
- **Δ > 0** : Régime **immiscible** (séparation de phases)

### Piège de confinement
```
V(x, t) = ½ m ω_final² x² + V₀(t) exp(−x²/2w²)
```

où **V₀(t)** décroît progressivement selon une rampe cosinus :
```
V₀(t) = V_barrière × ½ [1 + cos(π τ)]  avec τ ∈ [0, 1]
```

---

## 🛠️ Installation

### Prérequis

- Python ≥ 3.8
- NumPy ≥ 1.20
- Matplotlib ≥ 3.3
- SciPy (optionnel, pour analyses avancées)

---

## 🚀 Utilisation

### Exécution de la simulation
```bash
python bpm.py test_SI 1D
```

**Explication** :
- `bpm.py` : moteur de simulation (split-step BPM)
- `test_SI` : fichier de paramètres (sans extension `.py`)
- `1D` : module de grille et FFT (fichier `1D.py`)

### Modifier les paramètres

Éditez le fichier `test_SI.py` pour changer :

#### Espèces atomiques (ligne 57-68)
```python
# Cas immiscible : Rb-87 + Cs-133
g11_base = g1D_Rb_SI   / (hbar * omega_z_SI * L_scale_Rb)
g22_base = g1D_Cs_SI   / (hbar * omega_z_SI * L_scale_Cs)
g12_base = g1D_RbCs_SI / (hbar * omega_z_SI * 0.5 * (L_scale_Rb + L_scale_Cs))

# Cas miscible : Na-23 + Rb-87
g11_base = g1D_Na_SI   / (hbar * omega_z_SI * L_scale_Na)
g22_base = g1D_Rb_SI   / (hbar * omega_z_SI * L_scale_Rb)
g12_base = g1D_NaRb_SI / (hbar * omega_z_SI * 0.5 * (L_scale_Na + L_scale_Rb))
```

#### Facteur d'amplification (ligne 122)
```python
factor_interaction = 15.0  # Augmente la visibilité des interactions
```

#### Durée et résolution (ligne 135-145)
```python
tmax_SI = 1.500   # Durée totale (secondes)
dt_SI = 0.05e-3   # Pas de temps (50 µs)
images = 400      # Nombre d'images sauvegardées
```

---

## 📁 Structure du code
```
.
├── bpm.py                      # Moteur principal (split-step BPM)
├── 1D.py                        # Module grille, FFT, visualisation
├── test_SI.py                   # Paramètres physiques (Rb-Cs, Na-Rb, etc.)
├── README.md                    # Ce fichier
├── examples1D/                  # Dossier de sortie
│   └── test_SI/
│       ├── sectx_total.png      # Densité totale |ψ₁|² + |ψ₂|²
│       ├── sectx_psi1.png       # Densité Rb-87
│       ├── sectx_psi2.png       # Densité Cs-133
│       └── final_psi.npy        # État final (NumPy array)
└── docs/                        # Documentation (optionnel)
```

### Fichiers principaux

| Fichier | Description |
|---------|-------------|
| `bpm.py` | Boucle de propagation temporelle, appels FFT |
| `1D.py` | Grille spatiale, opérateurs, fonctions d'affichage |
| `test_SI.py` | Constantes physiques, couplages, conditions initiales |

---

## 🔬 Paramètres physiques

### Longueurs de diffusion (valeurs expérimentales)

| Interaction | aₛ (a₀) | Source |
|-------------|---------|--------|
| **Rb-Rb** | 94.7 | [van Kempen et al., PRL 2002] |
| **Cs-Cs** | 280 | [Chin et al., RMP 2010] |
| **Rb-Cs** | 650 | [Feshbach resonance, PRL 2005] |
| **Na-Na** | 52.9 | [Crubellier et al., EPJD 1999] |
| **Na-Rb** | 30 | [Estimation théorique] |

### Piège harmonique

- **Fréquence axiale** : ω_z = 2π × 6.8 Hz (faible confinement)
- **Fréquence transverse** : ω_⊥ = 2π × 160 Hz (fort confinement → quasi-1D)

### Unités de normalisation

- **Longueur** : a_z(Rb) = 4.14 µm
- **Temps** : τ₀ = 1/ω_z = 23.4 ms
- **Énergie** : E₀ = ℏω_z ≈ 4.5 × 10⁻³³ J

---

## 📊 Résultats

### Graphiques générés

Après exécution, 3 fichiers PNG sont créés dans `examples1D/test_SI/` :

1. **`sectx_total.png`** : Évolution spatio-temporelle de la densité totale |ψ₁|² + |ψ₂|²
2. **`sectx_psi1.png`** : Densité de l'espèce 1 (Rb-87 ou Na-23), colormap bleue
3. **`sectx_psi2.png`** : Densité de l'espèce 2 (Cs-133 ou Rb-87), colormap verte

**Normalisation** : Les 3 graphiques utilisent la même échelle [0, 1] pour comparaison directe.

### Cas typiques observés

| Système | Δ | Comportement |
|---------|---|--------------|
| **Rb-87 + Cs-133** | +39.5 | Immiscible : interface nette, pas de traversée |
| **Na-23 + Rb-87** | −0.8 | Miscible : collision quasi-élastique |
| **Li-7 + Na-23** | N/A | Instable (aₛ(Li) < 0 → collapse) |

---

## 📚 Références

### Articles scientifiques

1. **Dalibard, J.** (2024). *Solitons et ondes de matière*, Collège de France  
   [Lien cours](https://www.college-de-france.fr/)

2. **Olshanii, M.** (1998). *Atomic scattering in the presence of an external confinement and a gas of impenetrable bosons*, PRL **81**, 938  
   [DOI: 10.1103/PhysRevLett.81.938](https://doi.org/10.1103/PhysRevLett.81.938)

3. **Ao, P. & Chui, S.T.** (1998). *Binary Bose-Einstein condensate mixtures in weakly and strongly segregated phases*, PRA **58**, 4836  
   [DOI: 10.1103/PhysRevA.58.4836](https://doi.org/10.1103/PhysRevA.58.4836)

4. **Schreck, F. et al.** (2002). *Formation of a matter-wave bright soliton*, Science **296**, 1290  
   [DOI: 10.1126/science.1071021](https://doi.org/10.1126/science.1071021)

5. **Busch, T. & Anglin, J.R.** (2000). *Motion of dark solitons in trapped Bose-Einstein condensates*, PRL **84**, 2298  
   [DOI: 10.1103/PhysRevLett.84.2298](https://doi.org/10.1103/PhysRevLett.84.2298)

### Ressources complémentaires

- **Pethick & Smith** (2002). *Bose-Einstein Condensation in Dilute Gases*, Cambridge University Press
- **Pitaevskii & Stringari** (2016). *Bose-Einstein Condensation and Superfluidity*, Oxford University Press

---

## 👥 Auteurs

**Wan Penet** & **Maxime Foucher**  
*Projet de Recherche & Développement - LYRIDS - ECE Lyon, 2025/2026*

**Superviseurs** :  
- **M. Laurent Delisle** (Professeur chercheur ECE Paris)  
- **M. Amine Jaouadi** (Professeur chercheur ECE Paris)

---

## 📧 Contact

Pour toute question ou suggestion :  
📧 **wan.penet@edu.ece.fr**  
📧 **maxime.foucher@edu.ece.fr / myximefoucher@outlook.fr**

---

## 📝 Licence

Ce projet est sous licence **MIT**. 

---

## 🙏 Remerciements

- **Jean Dalibard** (Collège de France) pour les notes de cours sur les BEC
- **Edgar Figueiras et al.** (Université de Vigo) pour le code de base BPM
- **Communauté NumPy/SciPy** pour les outils de calcul scientifique
- **Laurent DELISLE, Amine JAOUADI et Wan PENET**

---

### Applications potentielles

- **Qubits solitoniques** pour le calcul quantique tolérant aux fautes
- **Métrologie atomique** (interféromètres à ondes de matière)
- **Étude de la turbulence quantique** dans les superfluides

---

**⭐ Si ce projet vous est utile, n'hésitez pas à mettre une étoile sur GitHub !**
