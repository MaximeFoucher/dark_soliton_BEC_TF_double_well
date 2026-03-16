# File: test_potentiel_laurent_SI.py
# Run as: python bpm2.py test_dimensione 1D
#
# SYSTÈME D'UNITÉS SI PURES
# - Longueur: mètres (m)
# - Temps: secondes (s)
# - Énergie: joules (J)
# - Fonction d'onde: m^(-1/2) en 1D

import numpy as np

# ============================================================================
# 1. CONSTANTES PHYSIQUES ET PARAMÈTRES ATOMIQUES (SI)
# ============================================================================
hbar = 1.054571596e-34
uaumass = 1.66053873e-27
epsilon0 = 8.854187e-12

# Paramètres atomiques (Rubidium 87)
mass = 86.909 * uaumass
a_0_Bohr = 5.29177e-11
a_s = 94.7 * a_0_Bohr # Longueur de diffusion

# Fréquences du piège (SI)
omega_z_SI = 2 * np.pi * 6.8 # Fréquence axiale (faible)
omega_rho_SI = 2 * np.pi * 160.0 # Fréquence radiale (forte)

# ============================================================================
# 2. DÉFINITION DES ÉCHELLES DE RÉFÉRENCE (SCALING)
# ============================================================================
# on calcule les facteurs pour convertir SI -> Moteur
# L_scale : Unité de longueur (~ 4 microns)
# T_scale : Unité de temps (~ 23 ms)
# E_scale : Unité d'énergie

L_scale = np.sqrt(hbar / (mass * omega_z_SI))
T_scale = 1.0 / omega_z_SI
E_scale = hbar * omega_z_SI

print("="*70)
print(f"PARAMÈTRES PHYSIQUES & CONVERSION")
print(f" Unité Longueur (a_z) : {L_scale*1e6:.2f} µm")
print(f" Unité Temps (1/w_z) : {T_scale*1000:.2f} ms")
print("="*70)

# ============================================================================
# 3. PARAMÈTRES DE SIMULATION (Définis en SI -> Convertis)
# ============================================================================
Nx = 1024 # Nombre de points de grille (puissance de 2 préférée)
Ny = Nx

# --- ESPACE ---
# On veut une fenêtre de +/- 150 microns (réaliste pour un BEC)
xmax_SI = 150e-6
ymax_SI = xmax_SI
xmax = xmax_SI / L_scale
ymax = xmax # On donne au moteur la valeur convertie (~35)

# --- TEMPS ---
# Durée totale : 500 ms
tmax_SI = 0.500
tmax = tmax_SI / T_scale # Conversion (~ 21 unités de temps)

# Pas de temps : 0.05 ms
dt_SI = 0.05e-3
dt = dt_SI / T_scale # Conversion (~ 0.002)

images = 200
output_choice = 3
fixmaximum = 0
absorb_coeff = 0.0

# ============================================================================
# 4. INTERACTIONS ET PHYSIQUE (SI -> Convertis)
# ============================================================================
# Calcul du g1D physique réel (Joules * m)
g1D_SI = 2 * hbar * omega_rho_SI * a_s

# Conversion en "g" sans dimension pour le code
g_simu = g1D_SI / (E_scale * L_scale)

# On applique aux espèces
g11 = g_simu * 1.5 # Un peu plus fort pour bien voir la répulsion
g22 = g_simu * 1.5
g12 = g_simu * 0.8 # Miscible (0.8 < 1.0)
g21 = g12

# Potentiel Chimique cible (en unités hw)
mu = mu_Rb = mu_Cs =15.0

# --- PARAMÈTRES DU MOUVEMENT (SI) ---
d_well_SI = 60e-6 # Distance initiale des puits : 60 µm
v_merge_SI = 100e-6 # Vitesse de rapprochement : 100 µm/s

# Conversion pour le moteur
d_well = d_well_SI / L_scale
v_merge = v_merge_SI / (L_scale/T_scale)

print(f"SIMULATION CONFIGURÉE:")
print(f" xmax : +/- {xmax_SI*1e6:.0f} µm")
print(f" tmax : {tmax_SI*1e3:.0f} ms")
print(f" Distance Puits : +/- {d_well_SI*1e6:.0f} µm")
print("="*70)

# Paramètres ajustés pour la fusion
v0_barriere = 20.0    # Hauteur de la barrière au centre à t=0
w_barrier = 8.0       # Largeur de la séparation (ajustez pour écarter les puits)
omega_final = 0.2     # Définit la courbure du piège final (le "serrage")

t_init_descente = 0
t_descente = tmax/2

def get_V0(t):
    """Fait descendre la barrière vers 0."""
    if t < t_init_descente:
        return v0_barriere
    elif t > t_init_descente + t_descente:
        return 0.0
    else:
        # Descente fluide
        progress = (t - t_init_descente) / t_descente
        smooth = 0.5 * (1 - np.cos(np.pi * progress))
        return v0_barriere * (1 - smooth)

def V_trap(x, t):
    # 1. Le piège final (parabole fixe qui maintient le condensat serré)
    V_harmonic = 0.5 * (omega_final**2) * (x**2)
    
    # 2. La barrière séparatrice qui s'efface
    V_barrier = get_V0(t) * np.exp(-(x**2) / (2 * w_barrier**2))
    
    return V_harmonic + V_barrier


# ============================================================================
# 6. FONCTION D'ONDE INITIALE (t=0)
# ============================================================================
def psi_0(x, y):
# Paramètres des Solitons (Trous)
    width_barrier = 1.0 # Largeur ~ 4 µm
    depth_barrier = 1.0 # Profondeur max (vide total)

    # Position initiale des puits (t=0)
    # Les solitons sont placés au fond des puits
    x_pos = d_well
    # Calcul du potentiel initial
    V0 = V_trap(x, 0.0)
    # --- Densité Thomas-Fermi ---
    n_TF_left = np.maximum((mu - V0) / g11, 0)
    n_TF_right = np.maximum((mu - V0) / g22, 0)
    # --- Création des Trous (Solitons) ---
    barrier_left = 1.0 - depth_barrier * np.exp(-(x + x_pos)**2 / (2 * width_barrier**2))
    barrier_right = 1.0 - depth_barrier * np.exp(-(x - x_pos)**2 / (2 * width_barrier**2))
    # --- Assemblage ---
    psi1 = np.sqrt(n_TF_left) * barrier_left * (x < 0)
    psi2 = np.sqrt(n_TF_right) * barrier_right * (x > 0)
    return np.array([psi1, psi2], dtype=complex)


# ============================================================================
# 7. ÉVOLUTION TEMPORELLE
# ============================================================================
def V(x, y, t, psi):

    psi1 = psi[0]
    psi2 = psi[1]

    V1 = V_trap(x, t) + g11 * np.abs(psi1)**2 + g12 * np.abs(psi2)**2
    V2 = V_trap(x, t) + g22 * np.abs(psi2)**2 + g21 * np.abs(psi1)**2
    return np.array([V1, V2])