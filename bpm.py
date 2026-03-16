import numpy as np
import sys
import os
import importlib
import glob

# Preliminaries (handling directories and files)
sys.path.insert(0, './examples'+sys.argv[2])		  # adds to path the directory with examples
output_folder = './examples'+sys.argv[2]+'/'+sys.argv[1]  # directory for images and video output
if not os.path.exists(output_folder):                     # creates folder if it does not exist
    os.makedirs(output_folder)

try:              # Erase all image files (if exist) before starting computation and generating new output
    for filename in glob.glob(output_folder+'/*.png') :   
        os.remove( filename )
except: 
    pass

my = importlib.__import__(sys.argv[1])		     # imports the file with the details for the computation
build = importlib.__import__(sys.argv[2])	     # selects 1D or 2D

# Initialization of the computation

x, y = build.grid(my.Nx,my.Ny,my.xmax,my.ymax)		# builds spatial grid
psi = my.psi_0(x,y) 					# loads initial condition


L = build.L(my.Nx,my.Ny,my.xmax,my.ymax)		# Laplacian in Fourier space
linear_phase = np.fft.fftshift(np.exp(1.j*L*my.dt/2))            	# linear phase in Fourier space (including point swap)
border = build.absorb(x,y,my.xmax,my.ymax,my.dt,my.absorb_coeff)    # Absorbing shell at the border of the computational window

# CORRECTION: Dimensions correctes pour le tableau de sauvegarde
# (images+1) lignes (temps) × Nx colonnes (espace)
savepsi=np.zeros((my.images+1, my.Nx))     # Creates a matrix to save the data of |psi|^2 for the final plot
savepsi1=np.zeros((my.images+1, my.Nx))    # Creates a matrix to save the data of |psi1|^2 for the final plot (composante 1)
savepsi2=np.zeros((my.images+1, my.Nx))    # Creates a matrix to save the data of |psi2|^2 for the final plot (composante 2)
steps_image=int(my.tmax/my.dt/my.images)  # Number of computational steps between consecutive graphic outputs


# Main computational loop
print("calculating", end="", flush=True)
for j in range(steps_image*my.images+1):		# propagation loop
	if j%steps_image == 0:  # Generates image output 
		build.output(x,y,psi,int(j/steps_image),j*my.dt,output_folder,my.output_choice,my.fixmaximum)
		# CORRECTION: Stockage correct avec les bonnes dimensions
		savepsi[int(j/steps_image), :] = build.savepsi(my.Ny,psi)
		savepsi1[int(j/steps_image), :] = build.savepsi1(my.Ny,psi)
		savepsi2[int(j/steps_image), :] = build.savepsi2(my.Ny,psi)
		print(".", end="", flush=True)

	V = my.V(x,y,j*my.dt,psi)			# potential operator
	psi *= np.exp(-1.j*my.dt*V)			# potential phase
	
	if sys.argv[2] == "1D":
			# Système binaire: transformation de chaque composante
			psi[0] = np.fft.fft(psi[0])
			psi[1] = np.fft.fft(psi[1])
			psi[0] *= linear_phase
			psi[1] *= linear_phase
			psi[0] = border * np.fft.ifft(psi[0])
			psi[1] = border * np.fft.ifft(psi[1])
	else: 
		print("Not implemented")

# Final operations
# Generates some extra output after the computation is finished and save the final value of psi:

# ============================================================================
# CALCUL DU MAXIMUM GLOBAL pour normalisation cohérente entre les 3 graphiques
# ============================================================================
max_global = max(np.max(savepsi), np.max(savepsi1), np.max(savepsi2))
print(f"\nMaximum de densité global : {max_global:.4f}")

build.final_output(output_folder,x,steps_image*my.dt,psi,savepsi,my.output_choice,my.images,my.fixmaximum, True, True, max_global)
build.final_output(output_folder,x,steps_image*my.dt,psi,savepsi1,my.output_choice,my.images,my.fixmaximum, True, False, max_global)
build.final_output(output_folder,x,steps_image*my.dt,psi,savepsi2,my.output_choice,my.images,my.fixmaximum, False, True, max_global)
build.generate_plotly_plots(output_folder, x, my.images, steps_image*my.dt, savepsi)
print()