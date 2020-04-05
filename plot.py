## Python plotting routines for MCPhac

import numpy as np
import matplotlib.pyplot as plt

# Physical constants (cgs)
h = 6.626075e-27  # erg s
k = 1.381e-16     # erg K-1
c = 2.998e10      # cm s-1
me = 9.10939e-28  # g
keV = 1.602177e-9 # 1000 erg

# Default parameters
it_def = 6 # number of iterations
N_def = 200 # number of depth points=

# Import run parameters
def load_params():
    with open('OUT/runparams','r') as f:
        # Ndepths,Ndepthsnu,Nmu,Nfreq,Maxfactor,Maxfactornu,MaxfracTch,Maxcoltau,Tguess,Maxiter,anist = [0 for i in range(15)]        	
    # for (line,var) in zip(f,(Ndepths,Ndepthsnu,Nmu,Nfreq,Maxfactor,Maxfactornu,MaxfracTch,Maxcoltau,Tguess,Maxiter,anist)):
        params = {}
        for line in f:
            params[line.split()[0]] = eval(line.split()[1])

    return params

#print(params)

def PlotSpectrum(it=it_def,N=N_def):

    with open('OUT/SurfaceFluxes.%d.%d.dat'%(N,it)) as f:
        nu,F = [],[]
        next(f)
        for line in f:
            nu.append(eval(line.split()[0]))
            F.append(eval(line.split()[1]))

    fig,ax = plt.subplots(1,1)
    ax.loglog(nu,F)
    ax.set_xlabel(r'$\nu$ (Hz)')
    ax.set_ylabel(r'$F$ (erg s$^{-1}$ cm$^{-2}$)')
    plt.show()

PlotSpectrum()
