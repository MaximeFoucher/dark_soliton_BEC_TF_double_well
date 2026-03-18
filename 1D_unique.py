import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm
import shutil
import platform


def grid(Nx,Ny,xmax,ymax):
	x = np.linspace(-xmax, xmax-2*xmax/Nx, Nx)     # x variable
	y = 0                 # not used, but a value must be given
	return x,y;

# Builds the Laplacian in Fourier space

def L(Nx,Ny,xmax,ymax):
	kx = np.linspace(-Nx/4/xmax, Nx/4/xmax-1/2/xmax, Nx)     # x variable
	return (2*np.pi*1.j*kx)**2 

# Introduces an absorbing shell at the border of the computational window

def absorb(x,y,xmax,ymax,dt,absorb_coeff):
	wx = xmax/20
	return np.exp(-absorb_coeff*(2-np.tanh((x+xmax)/wx)+np.tanh((x-xmax)/wx))*dt);

# Saves the data of abs(psi)**2 at different values of t

def savepsi(Ny,psi):
	return abs(psi)**2

# Defines graphic output: |psi|^2 is depicted

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
    
    # --- DENSITÉ NORMALISÉE entre 0 et 1 ---
    density = np.abs(psi)**2
    local_max = np.max(density)
    if local_max > 0:
        density_norm = density / local_max
    else:
        density_norm = density
    plt.plot(x, density_norm, label=r'$|\psi|^2$ (normalisé)', linewidth=2, color='blue')

    # --- AFFICHAGE DU POTENTIEL scalé sur [0, 1] ---
    v_ext = my.V_trap(x)
    v_display = np.copy(v_ext)
    v_limit = my.mu * 1.5  # On coupe à 150% du potentiel chimique
    v_display[v_display > v_limit] = v_limit
    # Scalé sur [0, 1] comme la densité normalisée
    v_scaled = v_display / v_limit
    plt.plot(x, v_scaled, 'r--', alpha=0.6, label='Double Well (scalé)')

    plt.xlabel('$x$')
    plt.ylabel('Densité normalisée / Potentiel')
    plt.title(f'$t = {t:.2f}$')
    plt.legend(loc='upper right')
    
    if fixmaximum > 0:
        plt.ylim(0, fixmaximum)
    else:
        plt.ylim(0, 1.4)  # fixe à 1.4 car le max normalisé vaut toujours 1

    if (output_choice==2) or (output_choice==3):
        plt.savefig(f'{folder}/fig{num}.png', dpi=100, bbox_inches='tight')

    if (output_choice==1) or (output_choice==3):
        plt.show(block=False)
        fig.canvas.flush_events()

# Some operations after the computation is finished: save the final value of psi, generate videos and builds
# the final plot: a contour map of the y=0 cut as a function of x and t

def final_output(folder,x,Deltat,psi,savepsi,output_choice,images,fixmaximum):

	np.save(folder,psi)		# saves final wavefunction

	if (output_choice==2) or (output_choice==3):
		movie(folder)	                        # creates video

	# Now we make a plot of the evolution depicting the 1D cut at y=0
	tvec=np.linspace(0,Deltat*images,images+1)
	tt,xx=np.meshgrid(tvec,x)
	figtx = plt.figure("Evolution of |psi(x)|^2")
	plt.clf()
	figtx.set_size_inches(8,6)

	# --- NORMALISATION GLOBALE entre 0 et 1 ---
	global_max = savepsi.max()
	if global_max > 0:
		toplot = savepsi / global_max
	else:
		toplot = savepsi

	if fixmaximum > 0:
		toplot[toplot > fixmaximum] = fixmaximum

	levels = np.linspace(0, 1, 100)
	plt.contourf(xx, tt, toplot, levels=levels, cmap=cm.jet, antialiased=False)
	cbar = plt.colorbar(ticks=[0, 0.25, 0.5, 0.75, 1.0])
	cbar.ax.set_yticklabels(['0.0', '0.25', '0.5', '0.75', '1.0'])
	plt.xlabel('$x (µm)$')
	plt.ylabel('$t (ms)$')
	cbar.set_label('$|\psi|^2$', fontsize=14)

	figname = folder+'/sectx.png'
	plt.savefig(figname)
	plt.show()


# Generates video from the saved figures. This function is called by final_output
def movie(folder):
	folder.replace('.','')
	examplename=folder[13:]

	video_options='vbitrate=4320000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq'		
	
	if platform.system() == 'Windows':
		try:
			shutil.copyfile('mencoder.exe', folder+'/mencoder.exe')
			os.chdir(folder)
			command ='mencoder "mf://fig*.png" -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie_'+examplename+'.avi'
			os.system(command)
			try:
				os.remove('mencoder.exe')
			except:
				print("Could not delete mencoder.exe in example directory")
			os.chdir('../../')
		except:
			print("Error making movie with mencoder in windows")
	else:
		try:
			command1 ='mencoder "mf://'+folder+'/fig*.png" -mf fps=25 -o /dev/null -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:'+video_options
			command2 ='mencoder "mf://'+folder+'/fig*.png" -mf fps=25 -o ./'+folder+'/movie_'+examplename+'.avi -ovc lavc -lavcopts vcodec=mpeg4:vpass=2:'+video_options
			os.system(command1)
			os.system(command2)
		except:
			print("Error making movie with mencoder in Linux")


## delete temporary files:
try:
	os.remove('divx2pass.log')
except:
	pass
			    
try:
	shutil.rmtree('__pycache__')
except:
	pass
		
try:
	shutil.rmtree('examples1D/__pycache__/')
except:
	pass