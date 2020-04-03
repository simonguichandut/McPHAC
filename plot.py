## Python plotting routines for MCPhac

import numpy as np
import matplotlib.pyplot as plt

# Default parameters
it_def = 6 # number of iterations
N_def = 200 # number of depth points

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
