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

def Fnu_norm(x):    # Normalized (x=hnu/kT)
    T = 1
    nu = kB*T*x/h
    return nu*Fnu(nu,T)/sigma/T**4

def run_McPHAC(logTeff=6.5,gsurf=2.43e14,logymin=-5.0,logymax=2.0,Ndepths=200,maxfactor=1,
            Ndepthsnu=200,maxfactornu=1,Nmu=10,Nfreq=100,maxfracTch=0.0001,maxiter=20,
            anist=0,maxcoltau=80.0,Tguess=0.264837817):
    """Can change any parameter from defaults (from source code .McPHAC.bash).
       Assuming code has been compiled (make McPHAC has been run)"""
    
    subprocess.call(["./McPHAC"]+[str(x) for x in(logTeff,gsurf,logymin,logymax,Ndepths,maxfactor,Ndepthsnu,maxfactornu,Nmu,Nfreq,maxfracTch,maxiter,anist,maxcoltau,Tguess)])
    

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
    """Finds number of iterations done by McPHAC (= # of TempProfile.dat files)"""

    i = 1
    file = get_path(run)+'TempProfile.'+str(load_params(run)['Ndepths'])+'.'+str(i)+'.dat'
    while os.path.exists(file):
        i+=1
        file = file[:-5] + str(i) + '.dat'
    return i-1


#------------------------------ Data file read functions -------------------------------

def read_datfile(file,columns):
    """Reads given columns indices of .dat file and returns them as numpy arrays"""
    arr = [[] for i in range(len(columns))]
    with open(file,'r') as f:
        next(f) # assuming first line of file should be skipped
        for line in f:
            for i,c in enumerate(columns):
                arr[i].append(eval(line.split()[c]))
    return (np.array(x) for x in arr)
    

def read_spectrum(run='current'):
    """Reads the (nu,Fnu) values from the last iteration SurfaceFluxes.dat file"""
    params = load_params(run)
    Flux_file = get_path(run) + ('SurfaceFluxes.%d.%d.dat'%(params['Ndepths'], num_iterations(run)))
    return read_datfile(Flux_file,columns=(0,1))

def read_Tprofile(run='current'):
    """Reads the (y,T) values from the last iteration TempProfile.dat"""
    params = load_params(run)
    Temp_file = get_path(run) + ('TempProfile.%d.%d.dat'%(params['Ndepths'], num_iterations(run)))
    return read_datfile(Temp_file,columns=(0,1))
        
def read_taufile(run='current'):
    """Reads the (E,y(tau=1)) values from tau.dat"""
    tau_file = get_path(run) + 'tau.dat'
    return read_datfile(tau_file,columns=(0,1))

#-------------------------------- General plotting functions --------------------------------

def NormSpectrum(run='current',ax=None,show=True,color='k',ls='-',label=None):
    """Normalized spectrum plot : nu*Fnu/sigmaTeff^4 vs hnu/kT. Option to pass axis handle as argument"""

    params = load_params(run)
    nu,Fnu = read_spectrum(run)

    if ax==None:
        _,ax = plt.subplots(1,1)
    ax.loglog(h*nu/kB/params['Teff'], nu*Fnu/sigma/params['Teff']**4,color=color,ls=ls,lw=0.8,label=label)
    
    if show:
        plt.show()
    else:
        return ax


def TempProfile(run='current',ax=None,show=True,color='k',ls='-',label=None):
    """Temperature (normalized to Teff) vs column depth plot"""

    params = load_params(run)
    y,T = read_Tprofile(run)

    if ax==None:
        _,ax = plt.subplots(1,1)
    ax.loglog(y,T/params['Teff'],color=color,ls=ls,lw=0.8,label=label)

    if show:
        plt.show()
    else:
        return ax




#-------------------------------- Specific plotting functions --------------------------------

def Make_fig1(force_recalculate=False):
    """ Figure 1 from Guichandut 2020. Same idea as fig 1 from Haakonsen et al. 2012 """

    # Create figure
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'$h\nu/kT_{eff}$',fontsize=14)
    ax.set_ylabel(r'$\nu F_\nu/\sigma T_{eff}^4$',fontsize=14)
    ax.set_xlim([0.1,100])
    ax.set_ylim([0.001,1])
    ax.tick_params(which='both',direction='in')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    
    # Plot Planck function
    x = np.logspace(-1,2,200)
    ax.loglog(x,Fnu_norm(x),'k-',label='BB',lw=0.8)
    
    # Compile code
    subprocess.call(["make","McPHAC"])

    # Run code at different effective temperatures and save the runs
    for logTeff, color in zip( (6.5,6.0,5.5), ('b','r','g') ):

        run_name = ("run_logTeff_%.1f_g14_2.4"%logTeff)

        # Don't run McPHAC if run is already saved (unless specified to do it anyway)
        if os.path.exists("OUT/" + run_name) and (not force_recalculate):
            pass
        else:
            run_McPHAC(logTeff=logTeff,gsurf=2.4e14)
            subprocess.call(["./save_run", run_name])

        # Plot spectrum
        NormSpectrum(run=run_name,ax=ax,show=False,color=color,label=(r'$T_{eff}=10^{%.1f}$K'%logTeff))

    ax.legend(frameon=False,loc='best')
    fig.savefig('figures/fig1.pdf')

def Make_Haakonsenfig1():
    """ Recreates figure 1 from Haakonsen et al. 2012 """

    # Create figure
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'$h\nu/kT_{eff}$',fontsize=14)
    ax.set_ylabel(r'$\nu F_\nu/\sigma T_{eff}^4$',fontsize=14)
    ax.set_xlim([0.1,100])
    ax.set_ylim([0.001,1])
    ax.tick_params(which='both',direction='in')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    # Plot Planck function
    x = np.logspace(-1,2,200)
    ax.loglog(x,Fnu_norm(x),'k-',label='BB',lw=0.8)

    # Compile code
    subprocess.call(["make","McPHAC"])

    # Run code at different effective temperatures and plot spectrum
    # Special linestyles: https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/linestyles.html
    for logTeff, ls in zip( (6.5,6.2,5.9,5.6,5.3), ('--',(0,(5,10)),':','-.',(0,(3,5,1,5))) ):
        run_McPHAC(logTeff=logTeff)
        NormSpectrum(ax=ax,show=False,ls=ls,label=(r'$T_{eff}=10^{%.1f}$K'%logTeff))

    ax.legend(frameon=False,loc='best')
    fig.savefig('figures/haakonsenfig1.pdf')


def Make_Haakonsenfig2():

    # Create figure
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'$h\nu (keV)$',fontsize=14)
    ax.set_ylabel(r'$F_{aniso}/F_{iso}$',fontsize=14)
    ax.set_xlim([0.1,10])
    ax.set_ylim([0.995,1.035])
    ax.tick_params(which='both',direction='in')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')

    # Compile code
    subprocess.call(["make","McPHAC"])

    # Run code at different effective temperatures and calculate ratios of fluxes
    for logTeff, ls in zip( (6.5,6.2,5.9,5.6,5.3), ('-','--',(0,(5,10)),':','-.')):

        # Run isotropic
        run_McPHAC(logTeff=logTeff,anist=0)
        nu,Fnu_iso = read_spectrum()

        # Run anisotropic
        run_McPHAC(logTeff=logTeff,anist=1)
        nu,Fnu_aniso = read_spectrum()

        # Plot
        ax.semilogx(h*nu/keV, Fnu_aniso/Fnu_iso, color='k', ls=ls, lw=0.8, label=(r'$T_{eff}=10^{%.1f}$K'%logTeff))

    ax.legend(frameon=False,loc='best')
    fig.savefig('figures/haakonsenfig2.pdf')










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
elif f == 'TempProfile': TempProfile(run)
elif f == 'Make_fig1': Make_fig1()
elif f == 'Make_Haakonsenfig1': Make_Haakonsenfig1()
elif f == 'Make_Haakonsenfig2': Make_Haakonsenfig2()






