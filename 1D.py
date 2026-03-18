import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm
import shutil
import platform
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


try:
    plt.rcParams.update({
    "text.usetex": False, # Désactive le vrai LaTeX
    "font.family": "STIXGeneral", # Utilise la police STIX (ressemble à Times/LaTeX)
    "mathtext.fontset": "stix", # Les maths en style LaTeX
    "font.size": 14,
    "axes.labelsize": 16,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    })
except Exception as e:
    print(f"Attention: Impossible de configurer LaTeX ({e}). Passage en mode standard.")

def grid(Nx,Ny,xmax,ymax):
    x = np.linspace(-xmax, xmax-2*xmax/Nx, Nx)
    y = 0
    return x,y

def L(Nx,Ny,xmax,ymax):
    kx = np.linspace(-Nx/4/xmax, Nx/4/xmax-1/2/xmax, Nx)
    return (2*np.pi*1.j*kx)**2 

def absorb(x,y,xmax,ymax,dt,absorb_coeff):
    wx = xmax/20
    return np.exp(-absorb_coeff*(2-np.tanh((x+xmax)/wx)+np.tanh((x-xmax)/wx))*dt)

def savepsi(Ny,psi):
    """
    Calcule la densité à sauvegarder pour le graphique final.
    Additionne les densités des deux composantes si le système est binaire.
    """
    if psi.ndim == 2:
        # Somme des densités |psi1|^2 + |psi2|^2
        return np.abs(psi[0])**2 + np.abs(psi[1])**2
    return np.abs(psi)**2

def savepsi1(Ny,psi):
    """
    Calcule la densité de la première composante (psi1) à sauvegarder pour le graphique final.
    Utile pour les systèmes binaires.
    """
    if psi.ndim == 2:
        return np.abs(psi[0])**2
    return np.abs(psi)**2

def savepsi2(Ny,psi):
    """
    Calcule la densité de la deuxième composante (psi2) à sauvegarder pour le graphique final."""
    if psi.ndim == 2:
        return np.abs(psi[1])**2
    return np.zeros_like(psi)  # Retourne un tableau de zéros si le système n'est pas binaire

def output(x,y,psi,n,t,folder,output_choice,fixmaximum):
    import sys
    import importlib
    example_name = sys.argv[1]
    my = importlib.import_module(example_name)
    importlib.reload(my) 

    if (output_choice==2) or (output_choice==3):
        num = str(int(n)).zfill(3)

    fig = plt.figure("1D plot")
    plt.clf()
    
    # --- GESTION BINAIRE (PSI1 et PSI2) ---
    if psi.ndim == 2:
        dens1 = np.abs(psi[0])**2
        dens2 = np.abs(psi[1])**2
        plt.plot(x, dens1, label=r'$|\psi_1|^2$', linewidth=2, color='blue')
        plt.plot(x, dens2, label=r'$|\psi_2|^2$', linewidth=2, color='green')
        current_max_density = max(np.max(dens1), np.max(dens2))
    else:
        density = np.abs(psi)**2
        plt.plot(x, density, label=r'$|\psi|^2$', linewidth=2, color='blue')
        current_max_density = np.max(density)

    # --- AFFICHAGE DU POTENTIEL ---
    v_ext = my.V_trap(x,t) # à modifier selon le fichier
    v_max_plot = current_max_density if current_max_density > 0 else 1.0
    

    v_display = np.copy(v_ext)
    # v_limit = my.mu * 1.5
    v_limit = np.maximum(my.mu_Rb, my.mu_Cs)
    v_display[v_display > v_limit] = v_limit
    
    v_scaled = (v_display / v_limit) * v_max_plot
    plt.plot(x, v_scaled, 'r--', alpha=0.6, label='Potentiel (scalé)')

    plt.xlabel('$x (µm)$')
    plt.ylabel('Densité / Potentiel')
    plt.title(f'$t = {t:.2f} ms$')
    plt.legend(loc='upper right')
    
    if fixmaximum > 0:
        plt.ylim(0, fixmaximum)
    else:
        plt.ylim(0, v_max_plot * 1.4)

    if (output_choice==1) or (output_choice==3):
        plt.show(block=False)
        fig.canvas.flush_events()

def final_output(folder,x,Deltat,psi,savepsi_data,output_choice,images,fixmaximum, boolean1, boolean2, max_global):
    """
    Génère le graphique de l'évolution spatio-temporelle (Contour Plot).
    La densité est normalisée entre 0 et 1 en utilisant max_global pour cohérence entre graphiques.
    
    Paramètres:
    -----------
    max_global : float
        Maximum de densité sur TOUS les graphiques (total, psi1, psi2)
        pour avoir la même échelle de colorbar 0-1 partout
    """
    # Sauvegarde de l'état final (première fois seulement)
    if boolean1 and boolean2:
        np.save(folder + '/final_psi', psi)

    if (output_choice==2) or (output_choice==3):
        if boolean1 and boolean2:  # Créer le film une seule fois
            movie(folder)

    # Création des vecteurs de temps et d'espace
    tvec = np.linspace(0, Deltat * images, savepsi_data.shape[0])
    tt, xx = np.meshgrid(tvec, x, indexing='ij')
    
    # ========================================================================
    # NORMALISATION : densité entre 0 et 1 avec le MAXIMUM GLOBAL
    # ========================================================================
    if max_global > 0:
        toplot = savepsi_data / max_global  # Normalisation 0-1
    else:
        toplot = savepsi_data
    
    # Application de fixmaximum si demandé (sur données normalisées)
    if fixmaximum > 0:
        # fixmaximum est interprété comme une fraction (0-1)
        toplot[toplot > fixmaximum] = fixmaximum
    
    # ========================================================================
    # CRÉATION DU GRAPHIQUE
    # ========================================================================
    fig = plt.figure(figsize=(12, 7))
    plt.clf()
    
    # Contour plot avec échelle 0-1 FIXE pour tous les graphiques
    levels = np.linspace(0, 1, 100)
    cs = plt.contourf(xx, tt, toplot, levels=levels, cmap=cm.jet, antialiased=False, extend='neither')
    
    cbar = plt.colorbar(cs, ticks=[0, 0.25, 0.5, 0.75, 1.0])
    cbar.ax.set_yticklabels(['0.0', '0.25', '0.5', '0.75', '1.0'])
    
    plt.xlabel('$x$ (µm)', fontsize=16)
    plt.ylabel('$t$ (ms)', fontsize=16)
    
    # Labels et titres selon le type de densité
    if boolean1 and boolean2:
        cbar.set_label('|ψ|²', fontsize=14)
        # title_str = 'Densité totale $|\psi_1|^2 + |\psi_2|^2$'
        suffix = 'total'
    elif boolean1:
        cbar.set_label('|ψ₁|²', fontsize=14)
        # title_str = 'Rubidium-87 : $|\psi_1|^2$'
        suffix = 'psi1'
    else:
        cbar.set_label('|ψ₂|²', fontsize=14)
        # title_str = 'Césium-133 : $|\psi_2|^2$'
        suffix = 'psi2'
    title_str = ''
    plt.title(title_str, fontsize=18, weight='bold')
    
    # Grille pour meilleure lisibilité
    plt.grid(True, alpha=0.2, linestyle='--')
    
    # ========================================================================
    # SAUVEGARDE AVEC NOM DE FICHIER ADAPTÉ
    # ========================================================================
    figname = f'{folder}/sectx_{suffix}.png'
    plt.savefig(figname, dpi=200, bbox_inches='tight')
    print(f"✓ {figname} — |ψ|² [0, 1]")
    
    plt.show()

def movie(folder):
    folder_clean = folder.replace('./','')
    examplename = folder_clean.split('/')[-1]
    video_options = 'vbitrate=4320000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq'        
    
    if platform.system() == 'Windows':
        try:
            if os.path.exists('mencoder.exe'):
                shutil.copyfile('mencoder.exe', folder+'/mencoder.exe')
                os.chdir(folder)
                command = 'mencoder "mf://fig*.png" -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie_'+examplename+'.avi'
                os.system(command)
                os.chdir('../../')
        except Exception as e:
            print(f"Error making movie: {e}")
    else:
        try:
            command1 = f'mencoder "mf://{folder}/fig*.png" -mf fps=25 -o /dev/null -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:{video_options}'
            command2 = f'mencoder "mf://{folder}/fig*.png" -mf fps=25 -o ./{folder}/movie_{examplename}.avi -ovc lavc -lavcopts vcodec=mpeg4:vpass=2:{video_options}'
            os.system(command1)
            os.system(command2)
        except:
            print("Error making movie with mencoder")

def generate_plotly_plots(folder, x, images, Deltat, savepsi_data, savepsi1_data, savepsi2_data):
    """
    Génère un graphique 3D interactif avec Plotly basé sur les données accumulées.
    """
    print("\n" + "="*70)
    print(" GÉNÉRATION DU GRAPHIQUE INTERACTIF PLOTLY (3D)")
    print("="*70)
    
    # Création des vecteurs de temps et d'espace
    tvec = np.linspace(0, Deltat * images, savepsi_data.shape[0])
    
    # Création de la figure 3D
    fig = go.Figure(data=[go.Surface(
        x=x, 
        y=tvec, 
        z=savepsi_data,
        colorscale='Jet',
        colorbar=dict(title="|ψ|²", thickness=20)
    )])

    fig.update_layout(
        title='Évolution Spatio-Temporelle (BPM 1D)',
        scene=dict(
            xaxis_title='x (µm)',
            yaxis_title='t (ms)',
            zaxis_title='|ψ|²',
            aspectmode='manual',
            aspectratio=dict(x=2, y=1.5, z=1),
            camera=dict(eye=dict(x=1.5, y=-1.5, z=1.3))
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )

    # Sauvegarde dans le dossier de simulation
    output_html = f"{folder}/simulation_interactive_3d.html"
    fig.write_html(output_html)
    print(f"✓ Fichier HTML généré : {output_html}")

    # Création des graphiques 3D pour psi1 

    # Création des vecteurs de temps et d'espace
    tvec = np.linspace(0, Deltat * images, savepsi1_data.shape[0])
    
    # Création de la figure 3D
    fig = go.Figure(data=[go.Surface(
        x=x, 
        y=tvec, 
        z=savepsi1_data,
        colorscale='Jet',
        colorbar=dict(title="|ψ1|²", thickness=20)
    )])

    fig.update_layout(
        title='Évolution Spatio-Temporelle (BPM 1D)',
        scene=dict(
            xaxis_title='x (µm)',
            yaxis_title='t (ms)',
            zaxis_title='|ψ1|²',
            aspectmode='manual',
            aspectratio=dict(x=2, y=1.5, z=1),
            camera=dict(eye=dict(x=1.5, y=-1.5, z=1.3))
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )

    # Sauvegarde dans le dossier de simulation
    output_html = f"{folder}/simulation_interactive_3d_psi1.html"
    fig.write_html(output_html)
    print(f"✓ Fichier HTML généré : {output_html}")

    # Création des graphiques 3D pour psi2

     # Création des vecteurs de temps et d'espace
    tvec = np.linspace(0, Deltat * images, savepsi2_data.shape[0])
    
    # Création de la figure 3D
    fig = go.Figure(data=[go.Surface(
        x=x, 
        y=tvec, 
        z=savepsi2_data,
        colorscale='Jet',
        colorbar=dict(title="|ψ2|²", thickness=20)
    )])

    fig.update_layout(
        title='Évolution Spatio-Temporelle (BPM 1D)',
        scene=dict(
            xaxis_title='x (µm)',
            yaxis_title='t (ms)',
            zaxis_title='|ψ2|²',
            aspectmode='manual',
            aspectratio=dict(x=2, y=1.5, z=1),
            camera=dict(eye=dict(x=1.5, y=-1.5, z=1.3))
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )

    # Sauvegarde dans le dossier de simulation
    output_html = f"{folder}/simulation_interactive_3d_psi2.html"
    fig.write_html(output_html)
    print(f"✓ Fichier HTML généré : {output_html}")

    

# Nettoyage
try:
    if os.path.exists('divx2pass.log'): os.remove('divx2pass.log')
    if os.path.exists('__pycache__'): shutil.rmtree('__pycache__')
    if os.path.exists('examples1D/__pycache__/'): shutil.rmtree('examples1D/__pycache__/')
except:
    pass