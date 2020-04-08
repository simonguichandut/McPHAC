## Python plotting routines for MCPhac

import sys
import os
import subprocess
import numpy as np
from scipy.special import expn
from scipy.integrate import quad

#--------------------------------- Plot settings ---------------------------------
# from https://jwalton.info/Embed-Publication-Matplotlib-Latex/

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

width_pt = 240 # MNRAS column width
inches_per_pt = 1 / 72.27
width_in = width_pt * inches_per_pt
golden_ratio = (5**.5 - 1) / 2
height_in = width_in * golden_ratio
fig_dim = (width_in,height_in)

fonts = {
        "axes.labelsize": 10,
        "font.size": 9,
        "legend.fontsize": 6,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "mathtext.default": "regular"
}
mpl.rcParams.update(fonts)



#------------------------------ Constants and basic functions ------------------------------

# (cgs)
h = 6.626075e-27  
kB = 1.381e-16     
c = 2.998e10     
arad = 7.5657e-15
sigma = 0.25*arad*c
keV = 1.602177e-9 

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
    ax.loglog(h*nu/kB/params['Teff'], nu*Fnu/sigma/params['Teff']**4,color=color,ls=ls,lw=0.9,label=label)
    
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
    ax.loglog(y,T/params['Teff'],color=color,ls=ls,lw=0.9,label=label)

    if show:
        plt.show()
    else:
        return ax

#-------------------------------- Specific plotting functions --------------------------------

def Make_fig12(force_recalculate=False):
    """ Figures 1 & 2 from Guichandut 2020. """

    # Create figures
    fig1,ax1 = plt.subplots(1,1,figsize=fig_dim)
    ax1.set_xlabel(r'$h\nu/kT_{eff}$')
    ax1.set_ylabel(r'$\nu F_\nu/\sigma T_{eff}^4$')
    ax1.set_xlim([0.1,100])
    ax1.set_ylim([0.001,1])
    ax1.tick_params(which='both',direction='in')
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')

    fig2,ax2 = plt.subplots(1,1,figsize=fig_dim)
    ax2.set_xlabel(r'$y$ (g cm$^{-2}$)')
    ax2.set_ylabel(r'$T/T_{eff}$')
    ax2.set_xlim([1e-6,10])
    ax2.set_ylim([0.1,100])
    ax2.tick_params(which='both',direction='in')
    ax2.yaxis.set_ticks_position('both')

    ax2b = ax2.twinx()
    # ax2b.tick_params(colors='b')
    ax2b.set_ylabel(r'$E_{\tau_\nu=1}$ (keV)',color='b')
    ax2.tick_params(which='both',direction='in')
    ax2b.set_ylim([1e-3,10])
    
    # Plot Planck function
    x = np.logspace(-1,2,200)
    ax1.loglog(x,Fnu_norm(x),'k-',label='blackbody',lw=0.8)
    
    # Compile code
    subprocess.call(["make","McPHAC"])

    # Run code at different effective temperatures and save the runs
    for logTeff, ls in zip( (6.5,6.0,5.5), ('--','-.',':') ):

        run_name = ("run_logTeff_%.1f_g14_2.4"%logTeff)

        # Don't run McPHAC if run is already saved (unless specified to do it anyway)
        if os.path.exists("OUT/" + run_name) and (not force_recalculate):
            pass
        else:
            run_McPHAC(logTeff=logTeff,gsurf=2.4e14,logymin=-7)
            subprocess.call(["./save_run", run_name])

        # Plot spectrum
        NormSpectrum(run=run_name,ax=ax1,show=False,ls=ls,label=(r'$T_{eff}=10^{%.1f}$K'%logTeff))

        # Plot temperature profile
        TempProfile(run=run_name,ax=ax2,show=False,ls=ls,label=(r'$T_{eff}=10^{%.1f}$K'%logTeff))

        # Add tau=1 data
        E,ytau1 = read_taufile(run=run_name)
        ax2b.loglog(ytau1,E,color='b',ls=ls,lw=0.9)

    ax1.legend(frameon=False,loc='best')
    ax2.legend(frameon=False,loc=2)
    fig1.tight_layout()
    fig2.tight_layout()
    fig1.savefig('figures/fig1.pdf', format='pdf', bbox_inches='tight')
    fig2.savefig('figures/fig2.pdf', format='pdf', bbox_inches='tight')

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
    fig.savefig('figures/haakonsenfig1.pdf', format='pdf', bbox_inches='tight')


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
    fig.savefig('figures/haakonsenfig2.pdf', format='pdf', bbox_inches='tight')










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
elif f == 'Make_fig12': Make_fig12()
elif f == 'Make_Haakonsenfig1': Make_Haakonsenfig1()
elif f == 'Make_Haakonsenfig2': Make_Haakonsenfig2()






