# PHYS642 : Radiative Processes in Astrophysics
## Simon Guichandut
### Course website : https://www.physics.mcgill.ca/~cumming/teaching/642/

This is my submission for the computationnal project. I use McPHAC ([Haakonsen et al. 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...749...52H/abstract)) 
to compute neutron star atmosphere spectra with different parameters and conditions, and analyze the results in a paper (Guichandut2020.pdf).

This repository is a clone of the original McPHAC repository ([link](https://github.com/McPHAC/McPHAC)). No changes were made to the source code.

This is the original readme: 

> The McGill Planar Hydrogen Atmosphere Code (McPHAC)

> To compile, make any needed changes to 'Makefile', and then type 'make McPHAC'.

> McPHAC takes commonly changed physical and computational parameters as
input arguments, and for convenience these can be set in in the script
'McPHAC.bash', which then allows running McPHAC with "./McPHAC.bash".

> Each function is contained in a separate file named after that
function, with the exception of 'McPHAC.c', which contains main().

Plotting routines were added in `plot.py`. Explanations of the output files were added in `OUTinfo.md`.

> python3 plot.py Make_figures

![](/figures/fig1.png)

![](/figures/fig2.png)

![](/figures/fig3.png) 

![](/figures/fig4.png) 
