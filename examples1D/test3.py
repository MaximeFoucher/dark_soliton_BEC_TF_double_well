# File: ./examples1D/testperso.py
# Run as    python bpm2.py test3 1D
# Génération de solitons sombres dans DEUX condensats de Bose-Einstein
# Chaque condensat a sa propre barrière pour générer des solitons

import numpy as np

# Define parameters
hbar=1.054571596e-34
clight=299792458.0
echarge=1.602176462e-19
emass=9.10938188e-31
pmass=1.67262158e-27
uaumass=1.66053873e-27
epsilon0=1.0e7/(4*np.pi*clight*clight)
kBoltzmann=1.3806503e-23

a_0 = 4 * np.pi * epsilon0 * hbar**2 / echarge / echarge / emass


# Define parameters
massRb = 86.909  # Atoms mass Cs 132.905 , Rb 86.909 (united atomic unit of mass)
massCs = 132.905  # Atoms mass Cs 132.905 , Rb 86.909 (united atomic unit of mass)
#mass  = mass * uaumass

Ntot= 20e4
omega_z = 2*np.pi*6.8
#U0 = 0.5 * mass * omega_rho**2
#U1 = 0.5 * mass * omega_z**2
#Period = 2*np.pi/(np.sqrt((2*U1/mass)))
#print(' Period =', Period)

#print('omega_rho =', omega_rho)
#print('omega_z =', omega_z)

#a_z = np.sqrt(hbar/mass/omega_z) 
a_s = 94.7*a_0 

g3d1 = (4 * np.pi * hbar**2 *a_s) / (massRb*uaumass) # page 19 Jean Dalibard
g3d2 = (4 * np.pi * hbar**2 *a_s) / (massCs*uaumass) # page 19 Jean Dalibard

print("g3d1 =", g3d1, "g3d2 =", g3d2)

Nx = 500                   # Grid points
Ny = Nx
dt = 0.02                  # Evolution step
tmax = 80                  # End of propagation
xmax = 50                   # x-window size
ymax = xmax                 # y-window size
images = 800                # number of .png images
absorb_coeff = 0.0          # Pas d'absorption - le piège confine naturellement
output_choice = 3           # 1: screen only, 2: save only, 3: both
fixmaximum = 0              # Fixes a maximum scale of |psi|**2 for the plots. If 0, it does not fix it.

# PARAMèTRES DU DOUBLE PUITS
v0_barriere = 15.0        # Hauteur de la barrière centrale du double puits
d_well = 30.0              # Demi-distance entre les deux puits


g11 = 1.5   # Densité maximale : nâ‚€ = Î¼/g donc si g augmente, pic diminue
g22 = 1.5  
g12 = 2.0   
g21 = g12

# PARAMèTRES PHYSIQUES DU PIèGE ET DU CONDENSAT
omega_trap = 0.1           
mu = mu_Rb = mu_Cs = 10.0                 

# Position réelle des minima du double puits (calculé à  partir de dV/dx = 0)
# Pour V(x) = V0 * [(x/d)^4 - 2*(x/d)^2 + 1], les minima sont à  x = Â±d
x_min_left = -d_well 
x_min_right = d_well 

R_TF = np.sqrt(2 * mu) / omega_trap
print(f"Rayon de Thomas-Fermi: R_TF = {R_TF:.2f}")
# print(f"Positions des minima du double puits: x_left = {x_min_left:.2f}, x_right = {x_min_right:.2f}")
print(f"Fréquence d'oscillation des solitons: omega_s = omega/sqrt(2) = {omega_trap/np.sqrt(2):.4f}")
print(f"Période d'oscillation: T = {2*np.pi/(omega_trap/np.sqrt(2)):.1f}")

t_init_descente = 0
t_descente = tmax/2
v_final = -15.0

# Paramètres ajustés pour la fusion
v0_barriere = 15.0    # Hauteur de la barrière au centre à t=0
w_barrier = 8.0       # Largeur de la séparation (ajustez pour écarter les puits)
omega_final = 0.1     # Définit la courbure du piège final (le "serrage")

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
    V_center = get_V0(t)
    V_barrier = V_center * np.exp(-(x**2) / (2 * w_barrier**2))
    
    return V_harmonic + V_barrier

# def get_V0(t):
#     t_start = tmax / 4
#     if t < t_start:
#         return v0_barriere
#     else:
#         # Baisse progressive vers un piège harmonique pur
#         progression = (t - t_start) / (tmax - t_start)
#         return v0_barriere * np.exp(-3 * progression)

# def V_trap(x, t=0):
#     # Confinement parabolique fixe (évite l'étalement)
#     V_harmonic = 0.5 * (0.1**2) * (x**2)
#     # Barrière centrale Gaussienne qui s'efface
#     V0 = get_V0(t)
#     w_barrier = 4.0
#     V_barrier = V0 * np.exp(-(x**2) / (2 * w_barrier**2))
#     return V_harmonic + V_barrier


# //////////////////////////////////////////////////////////////////////////////
# def get_V0(t):
#     
#     if t < t_init_descente:
#         return v0_barriere
#     
#     elif t > t_init_descente + t_descente:
#         return v_final
#     
#     else:
#         # Descente douce avec fonction cosinus
#         progress = (t - t_init_descente) / t_descente
#         smooth = 0.5 * (1 - np.cos(np.pi * progress))
#         return v0_barriere + smooth * (v_final - v0_barriere)

# def V_trap(x, t):
#     
#     v0 = get_V0(t)
#     return v0_barriere*(x/ d_well)**4 - 2 * v0 * (x / d_well)**2 + v0/2 + v0_barriere/2

# ============================================================================

# FONCTION D'ONDE INITIALE avec barrières pour générer les solitons
def psi_0(x, y):

    width_barrier = 1.0         # Largeur de la barrière (sigma)
    depth_barrier = 0.8         # Profondeur relative (0 à  1)

    
    # Position des barrières = centre des puits (minima du potentiel)
    x0_barrier_left = x_min_left    # Position calculée du minimum gauche
    x0_barrier_right = x_min_right  # Position calculée du minimum droit
    
    print(f"Barrières placées en x = {x0_barrier_left:.2f} (gauche) et x = {x0_barrier_right:.2f} (droite)")
    
    # CRà‰ATION DU CONDENSAT GAUCHE (Ïˆâ‚)
    # Potentiel vu par le condensat gauche (isolé)
    V_left = V_trap(x, t=0)
    n_left = (mu - V_left) / g11
    n_left = np.maximum(n_left, 0)
    
    # Barrière dans le condensat gauche
    # barrier_left = 1.0 - depth_barrier * np.exp(-(x - x0_barrier_left)**2 / (2 * width_barrier**2))
    barrier_left = barrier_right = 0.0

 
    
    # Condensat gauche = densité à— barrière à— masque spatial (x < 0)
    # psi1 = np.sqrt(n_left) * barrier_left * (x < 0)
    psi1 = np.sqrt(n_left) * (x < 0)
 

    
    # CRà‰ATION DU CONDENSAT DROIT (Ïˆâ‚‚)
    # Potentiel vu par le condensat droit (isolé)
    V_right = V_trap(x, t=0)
    n_right = (mu - V_right) / g22
    n_right = np.maximum(n_right, 0)
    
    # Barrière dans le condensat droit
    # barrier_right = 1.0 - depth_barrier * np.exp(-(x - x0_barrier_right)**2 / (2 * width_barrier**2))
    
    # Condensat droit = densité à— barrière à— masque spatial (x > 0)
    # psi2 = np.sqrt(n_right) * barrier_right * (x > 0)
    psi2 = np.sqrt(n_right) * (x > 0)
    
    return np.array([psi1, psi2], dtype=complex)



def V(x, y, t, psi):

    psi1 = psi[0]
    psi2 = psi[1]
    
    # Vâ‚ = V_trap + gâ‚â‚|Ïˆâ‚|Â² + gâ‚â‚‚|Ïˆâ‚‚|Â²
    V1 = V_trap(x, t) + g11 * np.abs(psi1)**2 + g12 * np.abs(psi2)**2
    
    # Vâ‚‚ = V_trap + gâ‚‚â‚‚|Ïˆâ‚‚|Â² + gâ‚‚â‚|Ïˆâ‚|Â²
    V2 = V_trap(x, t) + g22 * np.abs(psi2)**2 + g21 * np.abs(psi1)**2
    
    return np.array([V1, V2])