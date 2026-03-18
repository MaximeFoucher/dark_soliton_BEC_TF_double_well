# File: ./examples1D/test_BEC_unique.py
# Run as    python bpm_unique.py test_BEC_unique 1D
# Génération de solitons sombres dans un condensat de Bose-Einstein
# confiné dans un piège harmonique (profil Thomas-Fermi).
# Les solitons oscillent naturellement à la fréquence ω/√2.

import numpy as np


Nx = 500                   # Grid points
Ny = Nx
dt = 0.01                  # Evolution step (réduit pour meilleure précision)
tmax = 60                  # End of propagation (temps plus long pour voir oscillations)
xmax = 50                   # x-window size
ymax = xmax                 # y-window size
images = 400                # number of .png images
absorb_coeff = 0.0          # Pas d'absorption - le piège confine naturellement
output_choice = 1        # If 1, it plots on the screen but does not save the images
							# If 2, it saves the images but does not plot on the screen
							# If 3, it saves the images and plots on the screen
fixmaximum = 0            # Fixes a maximum scale of |psi|**2 for the plots. If 0, it does not fix it.
def V_trap(x):
    return 0.5 * omega_trap**2 * (x)**2 
#à modifier pour voir le x et la puissance 2


# ============================================================================
# PARAMÈTRES PHYSIQUES DU PIÈGE ET DU CONDENSAT
# ============================================================================
# Fréquence du piège harmonique (en unités naturelles), il dit à quel point les atomes sont confinés
omega_trap = 0.1           
# Potentiel chimique (détermine le nombre d'atomes / densité du condensat car n0 = mu/g)
mu = 10.0                 
# Force des interactions (g dans l'équation de Gross-Pitaevskii)
g_interaction = 1.5        # Coefficient de non-linéarité répulsive

# Rayon de Thomas-Fermi: R_TF = sqrt(2*mu / (m*omega^2))
# Dans nos unités: R_TF = sqrt(2*mu) / omega
R_TF = np.sqrt(2 * mu) / omega_trap
print(f"Rayon de Thomas-Fermi: R_TF = {R_TF:.2f}")
print(f"Fréquence d'oscillation des solitons: omega_s = omega/sqrt(2) = {omega_trap/np.sqrt(2):.4f}")
print(f"Période d'oscillation: T = {2*np.pi/(omega_trap/np.sqrt(2)):.1f}")


# FONCTION D'ONDE INITIALE 

#essayer avec potentiel gaussien
#essayer avec double well

def psi_0(x, y):
    
    # Profil de densité Thomas-Fermi
    # n(x) = (mu - 0.5*omega^2*x^2) / g si positif
    n_TF = (mu - V_trap(x)) / g_interaction
    n_TF = np.maximum(n_TF, 0)  # La densité ne peut pas être négative
    
    # Densité maximale au centre
    n0 = mu / g_interaction
    
    # ========================================================================
    # BARRIÈRE GAUSSIENNE pour générer les solitons
    # ========================================================================
    # Paramètres de la barrière
    x0_barrier = 0.0        # Position de la barrière (au centre)
    width_barrier = 1.0     # Largeur de la barrière (sigma) simgma = h_barre/sqrt(2*m*mu)
    depth_barrier = 0.8     # Profondeur relative (0 à 1)
                            # fait varier le nombre de solitons générés
    
    # La barrière crée une déplétion locale dans le condensat
    barrier = 1.0 - depth_barrier * np.exp(-(x - x0_barrier)**2 / (2 * width_barrier**2))
    
    # Fonction d'onde = sqrt(densité) * perturbation
    psi = np.sqrt(n_TF) * barrier
    
    
    return psi.astype(complex)


# POTENTIEL - Piège harmonique + Non-linéarité de Gross-Pitaevskii

def V(x, y, t, psi):
        
    # Non-linéarité répulsive (interactions entre atomes)
    V_interaction = g_interaction * np.abs(psi)**2
    
    # Potentiel total
    V_total = V_trap(x) + V_interaction
    
    return V_total

