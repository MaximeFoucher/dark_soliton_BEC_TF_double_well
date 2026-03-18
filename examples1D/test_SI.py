# File: test_SI.py
# Run as: python bpm.py test_SI 1D
#
# COUPLAGES VÉRIDIQUES calculés depuis g3D (Jean Dalibard)
# via intégration transverse (formule Olshanii)
#
# Système Rb-87 + Cs-133 en régime quasi-1D

import numpy as np


# ============================================================================
# 1. CONSTANTES PHYSIQUES
# ============================================================================
hbar = 1.054571596e-34
uaumass = 1.66053873e-27
echarge = 1.602176462e-19
emass = 9.10938188e-31
epsilon0 = 1.0e7 / (4 * np.pi * 2.998e8**2)
a_0 = 4 * np.pi * epsilon0 * hbar**2 / echarge**2 / emass  # Rayon de Bohr

# ============================================================================
# 2. PARAMÈTRES ATOMIQUES
# ============================================================================
mass_Rb = 86.909   # Rubidium-87 (u.a.m.)
mass_Cs = 132.905  # Césium-133 (u.a.m.)
massRb = mass_Rb*uaumass   # Rubidium-87 (u.a.m.)
massCs = mass_Cs*uaumass  # Césium-133 (u.a.m.)
mass_Li = 6.94
mass_Na = 22.99
massLi = mass_Li*uaumass  # Lithium-7 (u.a.m.) 
massNa = mass_Na*uaumass  # Sodium-23 (u.a.m.) 
mass_K = 39.10
massK = mass_K*uaumass  # Potassium-39 (u.a.m)



# Longueurs de diffusion (valeurs expérimentales)
a_s_Rb   = 94.7  * a_0   # Rb-Rb : 94.7 a₀ page 22 Ultracold rubidium interactions as=100 a₀ et page 43 as=90a0
a_s_Cs   = 280.0 * a_0   # Cs-Cs : ~280 a₀ page 1 High Resolution Feshbach Spectroscopy of Cesium as=280+-10ao
a_s_RbCs = 650.0 * a_0   # Rb-Cs : ~650 a₀ (Feshbach, valeur forte)
a_s_Li =  -27.3 * a_0   # Li-Li Abraham et al., PRL 74, 1315 (1995) 
a_s_Na = 52.9 * a_0    # Na-Na (valeur approximative)
a_s_LiNa =  22.0 * a_0    # Li-Na 
a_s_K =  174.0 * a_0    # K-K à modifier puis essayer avec le Rb et Na
a_s_NaK =  30.0 * a_0    # Na-K à modifier surement
a_s_NaRb =  30.0 * a_0    # Na-Rb à modifier surement
a_s_RbRb = a_s_Rb  # Rb-Rb (pour vérification cohérence)

# ============================================================================
# 3. PIÈGE HARMONIQUE
# ============================================================================
omega_z_SI   = 2 * np.pi * 6.8    # Fréquence axiale (faible)
omega_rho_SI = 2 * np.pi * 160.0  # Fréquence transverse (forte → quasi-1D)

# Échelles de référence pour l'adimensionnement
L_scale_Rb = np.sqrt(hbar / (massRb ) / omega_z_SI)  # ≈ 4.14 µm
L_scale_Cs = np.sqrt(hbar / (massCs ) / omega_z_SI)  # ≈ 3.34 µm
L_scale_Li = np.sqrt(hbar / (massLi ) / omega_z_SI)  # ≈ 6.63 µm
L_scale_Na = np.sqrt(hbar / (massNa ) / omega_z_SI)  # ≈ 5.17 µm
L_scale_K = np.sqrt(hbar / (massK ) / omega_z_SI)    # ≈ 4.90 µm
T_scale = 1.0 / omega_z_SI                                     # ≈ 23.4 ms
E_scale = hbar * omega_z_SI

print("="*75)
print("CALCUL DES COUPLAGES VÉRIDIQUES (Rb-87 + Cs-133)")
print("="*75)
print(f"Longueur oscillateur (Rb) : a_z = {L_scale_Rb*1e6:.3f} µm")
print(f"Longueur oscillateur (Cs) : a_z = {L_scale_Cs*1e6:.3f} µm")
print(f"Unité de temps            : τ₀ = {T_scale*1e3:.2f} ms")
print(f"Unité d'énergie           : E₀ = {E_scale/echarge*1e9:.2f} neV")
print("="*75)

# ============================================================================
# 4. CALCUL DES COUPLAGES (méthode correcte : g3D → g1D → adimensionné)
# ============================================================================

# --- Étape 1 : Couplage 3D (formule Jean Dalibard p.19) ---
print("\n--- ÉTAPE 1 : Couplage 3D (Jean Dalibard) ---")
g3d_Rb   = (4 * np.pi * hbar**2 * a_s_Rb)   / (massRb )
g3d_Cs   = (4 * np.pi * hbar**2 * a_s_Cs)   / (massCs)
g3d_RbCs = (4 * np.pi * hbar**2 * a_s_RbCs) / (np.sqrt(massRb * massCs) )
g3d_Li   = (4 * np.pi * hbar**2 * a_s_Li)   / (massLi)
g3d_Na   = (4 * np.pi * hbar**2 * a_s_Na)   / (massNa)
g3d_LiNa = (4 * np.pi * hbar**2 * a_s_LiNa) / (np.sqrt(massLi * massNa) )
g3d_K   = (4 * np.pi * hbar**2 * a_s_K)   / (massK)
g3d_NaK = (4 * np.pi * hbar**2 * a_s_NaK)   / (np.sqrt(massNa * massK) )
g3d_NaRb = (4 * np.pi * hbar**2 * a_s_NaRb)   / (np.sqrt(massNa * massRb) )
g3d_RbRb = (4 * np.pi * hbar**2 * a_s_RbRb) / (massRb )  # pour vérification cohérence

print(f"g3D(Rb-Rb) = {g3d_Rb:.4e} J·m³")
print(f"g3D(Cs-Cs) = {g3d_Cs:.4e} J·m³")
print(f"g3D(Rb-Cs) = {g3d_RbCs:.4e} J·m³")
print(f"g3D(Li-Li) = {g3d_Li:.4e} J·m³")
print(f"g3D(Na-Na) = {g3d_Na:.4e} J·m³")
print(f"g3D(Li-Na) = {g3d_LiNa:.4e} J·m³")
print(f"g3D(K-K) = {g3d_K:.4e} J·m³")
print(f"g3D(Na-K) = {g3d_NaK:.4e} J·m³")
print(f"g3D(Na-Rb) = {g3d_NaRb:.4e} J·m³")

# --- Étape 2 : Intégration transverse (Olshanii, PRL 81, 1998) ---
print("\n--- ÉTAPE 2 : Réduction quasi-1D (Olshanii) ---")
print("Formule : g₁D = g3D / (2π a_⊥²) = 2 ℏ ω_⊥ a_s")

a_perp_Rb = np.sqrt(hbar / (massRb ) / omega_rho_SI)
a_perp_Cs = np.sqrt(hbar / (massCs) / omega_rho_SI)

# Méthode directe (équivalente) page 19 du cours de Jean Dalibard
g1D_Rb_SI   = 2 * hbar * omega_rho_SI * a_s_Rb 
g1D_Cs_SI   = 2 * hbar * omega_rho_SI * a_s_Cs
g1D_RbCs_SI = 2 * hbar * omega_rho_SI * a_s_RbCs
g1D_Li_SI   = 2 * hbar * omega_rho_SI * a_s_Li
g1D_Na_SI   = 2 * hbar * omega_rho_SI * a_s_Na
g1D_LiNa_SI = 2 * hbar * omega_rho_SI * a_s_LiNa
g1D_K_SI    = 2 * hbar * omega_rho_SI * a_s_K
g1D_NaK_SI  = 2 * hbar * omega_rho_SI * a_s_NaK
g1D_NaRb_SI  = 2 * hbar * omega_rho_SI * a_s_NaRb
g1D_RbRb_SI = 2 * hbar * omega_rho_SI * a_s_RbRb  # pour vérification cohérence

print(f"g₁D(Rb-Rb) = {g1D_Rb_SI:.4e} J·m")
print(f"g₁D(Cs-Cs) = {g1D_Cs_SI:.4e} J·m")
print(f"g₁D(Rb-Cs) = {g1D_RbCs_SI:.4e} J·m")
print(f"g₁D(Li-Li) = {g1D_Li_SI:.4e} J·m")
print(f"g₁D(Na-Na) = {g1D_Na_SI:.4e} J·m")
print(f"g₁D(Li-Na) = {g1D_LiNa_SI:.4e} J·m")
print(f"g₁D(K-K) = {g1D_K_SI:.4e} J·m")
print(f"g₁D(Na-K) = {g1D_NaK_SI:.4e} J·m")
print(f"g₁D(Na-Rb) = {g1D_NaRb_SI:.4e} J·m")

# Vérification cohérence
verif_Rb = g1D_Rb_SI / (g3d_Rb / (2 * np.pi * a_perp_Rb**2))
print(f"\nVérification : g₁D(Rb) / [g3D(Rb) / (2π a_⊥²)] = {verif_Rb:.3f} (doit être ≈ 1)")

# --- Étape 3 : Adimensionnement ---
print("\n--- ÉTAPE 3 : Adimensionnement pour le code BPM ---")
# --------------------------------------- fonctinonne non miscible
# g11_base = g1D_Rb_SI   / (hbar * omega_z_SI * L_scale_Rb)
# g22_base = g1D_Cs_SI   / (hbar * omega_z_SI * L_scale_Cs)
# g12_base = g1D_RbCs_SI / (hbar * omega_z_SI * 0.5 * (L_scale_Rb + L_scale_Cs))
# g21_base = g12_base
# ---------------------------------------- foncionne miscible
g11_base = g1D_Na_SI   / (hbar * omega_z_SI * L_scale_Na)
g22_base = g1D_Rb_SI   / (hbar * omega_z_SI * L_scale_Rb)
g12_base = g1D_NaRb_SI / (hbar * omega_z_SI * 0.5 * (L_scale_Na + L_scale_Rb))
g21_base = g12_base
# ---------------------------------------- fonctionne pas
# g11_base = g1D_Li_SI   / (hbar * omega_z_SI * L_scale_Li)
# g22_base = g1D_Na_SI   / (hbar * omega_z_SI * L_scale_Na)
# g12_base = g1D_LiNa_SI / (hbar * omega_z_SI * 0.5 * (L_scale_Li + L_scale_Na))
# g21_base = g12_base
# ---------------------------------------- fonctionne pas 
# g11_base = g1D_K_SI   / (hbar * omega_z_SI * L_scale_K)
# g22_base = g1D_Na_SI   / (hbar * omega_z_SI * L_scale_Na)
# g12_base = g1D_NaK_SI / (hbar * omega_z_SI * 0.5 * (L_scale_Na + L_scale_K))
# g21_base = g12_base
# ----------------------------- test
# g11_base = g1D_Li_SI   / (hbar * omega_z_SI * L_scale_Li)
# g22_base = g1D_Li_SI   / (hbar * omega_z_SI * L_scale_Li)
# g12_base = g1D_Li_SI / (hbar * omega_z_SI * L_scale_Li)
# g21_base = g12_base


print(f"g₁₁ (base, Rb-Rb) = {g11_base:.6f}")
print(f"g₂₂ (base, Cs-Cs) = {g22_base:.6f}")
print(f"g₁₂ (base, Rb-Cs) = {g12_base:.6f}")

Delta_base = g12_base * g21_base - g11_base * g22_base
print(f"\nDiscriminant Δ (base) = {Delta_base:.6f}")
if Delta_base > 0:
    print(f"  → IMMISCIBLE")
else:
    print(f"  → Miscible")



Ntot = 1.2e4
# mu = 15
aohRb = np.sqrt(hbar/(massRb*omega_z_SI))
aohCs = np.sqrt(hbar/(massCs*omega_z_SI))
aohNa = np.sqrt(hbar/(massNa*omega_z_SI))
aohK = np.sqrt(hbar/(massK*omega_z_SI))
aohLi = np.sqrt(hbar/(massLi*omega_z_SI))


mu_temp_Rb = 0.5 * ((15*Ntot*a_s_Rb)/aohRb)**((2/5))* hbar*omega_z_SI  # estimation grossière pour le potentiel chimique à partir de la densité maximale
mu_temp_Cs = 0.5 * ((15*Ntot*a_s_Cs)/aohCs)**((2/5))*hbar*omega_z_SI  
mu_temp_Na = 0.5 * ((15*Ntot*a_s_Na)/aohNa)**((2/5))*hbar*omega_z_SI 
mu_temp_K = 0.5 * ((15*Ntot*a_s_K)/aohK)**((2/5))*hbar*omega_z_SI  
mu_temp_Li = 0.5 * ((15*Ntot*a_s_Li)/aohLi)**((2/5))*hbar*omega_z_SI 

mu_Rb = mu_temp_Rb/(hbar*omega_z_SI)  # estimation grossière pour le potentiel chimique à partir de la densité maximale
mu_Cs = mu_temp_Cs/(hbar*omega_z_SI)  
mu_Na = mu_temp_Na/(hbar*omega_z_SI)  
mu_K = mu_temp_K/(hbar*omega_z_SI)  
mu_Li = mu_temp_Li/(hbar*omega_z_SI) 



print(f"mu_Rb : ", mu_Rb, f"mu_Cs : ", mu_Cs)

# ============================================================================
# 5. FACTORISATION POUR VISIBILITÉ (optionnel)
# ============================================================================
print("\n" + "="*75)
print("STRATÉGIE : Factorisation pour rendre les interactions VISIBLES")
print("="*75)

# Choix : facteur 15 (compromis entre réalisme et visibilité)
factor_interaction = 15.0

g11 = factor_interaction * g11_base  # ≈ 0.86
g22 = factor_interaction * g22_base  # ≈ 3.13
g12 = factor_interaction* g12_base  # ≈ 6.49
g21 = g12

Delta = - g12 * g21 + g11 * g22
print(f"\nFacteur appliqué : {factor_interaction:.1f}×")
print(f"  (équivaut à ω_⊥ ~ 2π × {omega_rho_SI/(2*np.pi) * factor_interaction:.0f} Hz)")
print(f"\nCouplages effectifs dans le code :")
print(f"  g₁₁ (Rb-Rb) = {g11:.4f}")
print(f"  g₂₂ (Cs-Cs) = {g22:.4f}")
print(f"  g₁₂ (Rb-Cs) = {g12:.4f}")
print(f"\nDiscriminant Δ = {Delta:.4f}")
if Delta < 0:
    print(f"  → IMMISCIBLE (séparation visible)")
else:
    print(f"  → Miscible")
print("="*75)

# ============================================================================
# 6. PARAMÈTRES DE SIMULATION
# ============================================================================
Nx = 1024
Ny = Nx

# Espacea
xmax_SI = 150e-6  # ±150 µm
ymax_SI = xmax_SI
xmax = xmax_SI / L_scale_Rb
ymax = xmax

# Temps
tmax_SI = 1.500   # 500 ms
dt_SI = 0.05e-3   # 0.05 ms
tmax = tmax_SI / T_scale
dt = dt_SI / T_scale

images = 400
output_choice = 3
fixmaximum = 0
absorb_coeff = 0.0



# Paramètres double puits
d_well_SI = 60e-6
v_merge_SI = 100e-6
d_well = d_well_SI / L_scale_Rb
v_merge = v_merge_SI / (L_scale_Rb / T_scale)

print(f"\nSIMULATION :")
print(f"  Fenêtre spatiale : ±{xmax_SI*1e6:.0f} µm")
print(f"  Durée totale     : {tmax_SI*1e3:.0f} ms")
print(f"  Distance puits   : ±{d_well_SI*1e6:.0f} µm")
print("="*75)

# ============================================================================
# 7. POTENTIEL EXTERNE
# ============================================================================
v0_barriere = 5.0
w_barrier = 4.0
omega_final = 0.15

t_init_descente = 0
t_descente = tmax / 2



def get_V0(t):
    """Barrière centrale qui descend progressivement."""
    if t < t_init_descente:
        return v0_barriere
    elif t > t_init_descente + t_descente:
        return 0.0
    else:
        progress = (t - t_init_descente) / t_descente
        smooth = 0.5 * (1 - np.cos(np.pi * progress))
        return v0_barriere * (1 - smooth)

def V_trap(x, t):
    """
    Piège harmonique + barrière gaussienne centrale.
    """
    V_harmonic = 0.5 * (omega_final**2) * (x**2)
    V_barrier = get_V0(t) * np.exp(-(x**2) / (2 * w_barrier**2))
    return V_harmonic + V_barrier

# ============================================================================
# 8. CONDITION INITIALE (Thomas-Fermi)
# ============================================================================
def psi_0(x, y):
    """
    Deux condensats Thomas-Fermi séparés :
    - ψ₁ (Rb) à gauche (x < 0)
    - ψ₂ (Cs) à droite (x > 0)
    """
    V0 = V_trap(x, 0.0)
    
    # Densités Thomas-Fermi
    # ----------------------------------- fonctionne non miscible
    # n_TF_a = np.maximum((mu_Rb - V0) / g11, 0)
    # n_TF_b = np.maximum((mu_Cs - V0) / g22, 0)
    # ----------------------------------- fonctionne miscible
    n_TF_a = np.maximum((mu_Na - V0) / g11, 0)
    n_TF_b = np.maximum((mu_Rb - V0) / g22, 0)
    # ----------------------------------- deuxième cas
    # n_TF_a = np.maximum((mu_Li - V0) / g11, 0)
    # n_TF_b = np.maximum((mu_Na - V0) / g22, 0)
    # ----------------------------------- troisième cas
    # n_TF_a = np.maximum((mu_Na - V0) / g22, 0)
    # n_TF_b = np.maximum((mu_K - V0) / g11, 0)
    # -------------- test
    # n_TF_a = np.maximum((mu_Li - V0) / g11, 0)
    # n_TF_b = np.maximum((mu_Li - V0) / g22, 0)

    # Séparation spatiale initiale
    psi1 = np.sqrt(n_TF_a) * (x < 0)
    psi2 = np.sqrt(n_TF_b) * (x > 0)

    print(f"psi1 max (Rb) : {np.max(psi1):.4f}, psi2 max (Cs) : {np.max(psi2):.4f}")
    
    
    return np.array([psi1, psi2], dtype=complex)

# ============================================================================
# 9. POTENTIEL NON-LINÉAIRE (GPE)
# ============================================================================
def V(x, y, t, psi):
    """
    Potentiel effectif vu par chaque composante.
    """
    psi1, psi2 = psi[0], psi[1]
    
    V1 = V_trap(x, t) + g11 * np.abs(psi1)**2 + g12 * np.abs(psi2)**2
    V2 = V_trap(x, t) + g22 * np.abs(psi2)**2 + g21 * np.abs(psi1)**2
    
    return np.array([V1, V2])