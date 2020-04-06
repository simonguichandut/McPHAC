## Python plotting routines for MCPhac

import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expn
from scipy.integrate import quad

# Physical constants (cgs)
h = 6.626075e-27  # erg s
kB = 1.381e-16     # erg K-1
c = 2.998e10      # cm s-1
me = 9.10939e-28  # g
arad = 7.5657e-15
sigma = 0.25*arad*c
keV = 1.602177e-9 # 1000 erg


def Bnu(nu,T):    # Planck function
    return 2*h*nu**3/c**2 * 1/(np.exp(h*nu/kB/T) - 1)
def Fnu(nu,T):    # Flux coming out of isothermal atmosphere
    F = []
    for nui in nu:
        integrand = lambda t: Bnu(nui,T)*expn(2,t)
        F.append(2*np.pi*quad(integrand,0,np.inf)[0])
    return F

def Fnu_norm(x):    # Normalized
    T = 1
    nu = kB*T*x/h
    return nu*Fnu(nu,T)/sigma/T**4


#-------------------------------- Run info functions --------------------------------

def get_path(run='current'):
    """Path of saved files of given run"""
    path = 'OUT/' if run=='current' else ('OUT/'+run+'/') 
    return path

def load_params(run='current'):
    """Loads run paramaters saved in runparams file"""
    
    with open(get_path(run)+'runparams','r') as f:
        params = {}
        for line in f:
            params[line.split()[0]] = eval(line.split()[1])
    return params

def num_iterations(run='current'):
    """Finds number of iterations done by McPHAC (= # of JT.dat files)"""

    i = 1
    file = get_path(run)+'JT.'+str(i)+'.dat'
    while os.path.exists(file):
        i+=1
        file = get_path(run)+'JT.'+str(i)+'.dat'
    return i-1


#-------------------------------- General plotting functions --------------------------------

def NormSpectrum(run='current',ax=None,show=True,color='k',ls='-',label=None):
    """Normalized spectrum plot : nu*Fnu/sigmaTeff^4 vs hnu/kT. Option to pass axis handle as argument"""

    params = load_params(run)
    Flux_file = get_path(run) + ('SurfaceFluxes.%d.%d.dat'%(params['Ndepths'], num_iterations(run)))
    with open(Flux_file,'r') as f:
        nu,F = [],[]
        next(f)
        for line in f:
            nu.append(eval(line.split()[0]))
            F.append(eval(line.split()[1]))
        nu,F = np.array(nu),np.array(F)

    if ax==None:
        _,ax = plt.subplots(1,1)
    ax.loglog(h*nu/kB/params['Teff'], nu*F/sigma/params['Teff']**4,color=color,ls=ls,lw=0.8,label=label)
    
    if show:
        plt.show()
    else:
        return ax

def Otherplot():
    
    pass




#-------------------------------- Specific plotting functions --------------------------------

def Make_paperfig1():
    """ Recreates figure 1 from Haakonsen et al. 2012 """

    # Create figure
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'$h\nu/kT_{eff}$',fontsize=14)
    ax.set_ylabel(r'$\nu F_\nu/\sigma T_{eff}^4$',fontsize=14)
    ax.set_xlim([0.1,100])
    # ax.set_ylim([0.001,1])
    ax.set_ylim([0.01,1])
    ax.tick_params(which='both',direction='in')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    # Plot Planck function
    x = np.logspace(-1,2,200)
    ax.loglog(x,Fnu_norm(x),'k-',label='BB')

    # Compile code
    subprocess.call(["make","McPHAC"])

    # Default parameter values (from source code .McPHAC.bash)
    defaults = [6.5, 2.43e14, -5.0, 2.0, 200, 1, 200, 1, 10, 100, 0.0001, 20, 0, 80.0, 0.264837817]

    # Run code at different effective temperature and plot spectrum
    # Special linestyles: https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/linestyles.html
    for logTeff, ls in zip( (6.5,6.2,5.9,5.6,5.3), ('--',(0,(5,10)),':','-.',(0,(3,5,1,5))) ):
        cmnd = ["./McPHAC", str(logTeff)]
        cmnd += [str(x) for x in defaults[1:]]
        subprocess.call(cmnd)
        NormSpectrum(ax=ax,show=False,ls=ls,label=(r'$T_{eff}=10^{%.1f}K$'%logTeff))

    ax.legend(frameon=False,loc='best')
    plt.show()










## Command-line call
# Argument parsing : first arg (mandatory) is plotting function to run, second arg (optionnal) name of saved run
msg = "Command line call : $python3 plot.py [function_name] [run_name]"
if len(sys.argv) == 1 or len(sys.argv)>3:
    sys.exit(msg)
elif len(sys.argv) == 2:
    run = 'current' # No run name provided, use current run (saved by default in OUT/)
else:
    run = sys.argv[2]

f = sys.argv[1]
if f == 'NormSpectrum': NormSpectrum(run)
elif f == 'Make_paperfig1': Make_paperfig1()







