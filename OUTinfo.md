# Description of output files (OUT/)
Referencing equations from [Haakonsen et al. 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...749...52H/abstract)
 
# Main Results

### **[EmergentSpectrum.X.Y.dat](./OUT/EmergentSpectrum.200.1.dat)**
Spectrum of the emergent radiation, with angle dependence

* X : number of depth points
* Y : iteration number
* Columns:

    |  |  |
    | --- | --- |
    1) nu | Frequency (Hz)
    2) mu | Angle cosine
    3) I | Surface specific intensity at mu and nu

* Blocks : One per frequency, Nmu lines each
* Source code : [OutputSpectrum.c](./OutputSpectrum.c)

### **[SurfaceFluxes.X.Y.dat](./OUT/SurfaceFluxes.200.1.dat)**
Total emergent flux as a function of frequency

* X : number of depth points
* Y : iteration number
* Columns:

    |  |  |
    | --- | --- |
    1) nu | Frequency (Hz)
    2) F | Surface flux (erg/cm2/s)

* Source code : [OutputFluxes.c](./OutputFluxes.c)
* [SurfaceFluxest.X.Y.dat](./OUT/SurfaceFluxest.200.1.dat) : same structure, after temperature correction prodcedure

### **[TempProfile.X.Y.dat](./OUT/TempProfile200.1.dat)**
Temperature profile of the atmosphere

* X : number of depth points
* Y : iteration number
* Columns:
        
    |  |  |
    | --- | --- |
    1) y | Column depth (g/cm2)
    2) T | Temperature (K) (10^logT)
    3) T | Temperature (K)
    3) Tlogint | ?
    
* Source code : [PrintInterpolT.c](./PrintInterpolT.c)


## Other


### **[Temperatures.X.Y.dat](./OUT/Temperatures.200.1.dat)**

* X : number of depth points
* Y : iteration number
* Columns:
        
    |  |  |
    | --- | --- |
    1) #d| X
    2) logy | Column depth
    3) oldT | Old temperature 
    4) newT | New temperature
    5) deltaT | Difference
    6) deltaT/T | 
    7) avgf | ?
    8) fabs(avgf/sigmaf) | ?
    
* Blocks : One per frequency, X lines each
* Source code : [DeltaT.c](./DeltaT.c)


### **[Fluxes.X.Y.dat](./OUT/Fluxes.200.1.dat)**

* X : number of depth points
* Y : iteration number
* Columns:

    |  |  |
    | --- | --- |
    1) y | Column depth (g/cm2)
    2) T | Temperature (K)
    3) Teff | Effective temperature (constant)
    (weird name for actual contents?)

* Source code : [OutputFluxes.c](./OutputFluxes.c)
* [Fluxest.X.Y.dat](./OUT/Fluxest.200.1.dat) : same structure, after temperature correction prodcedure



### **[SurfaceFluxesNorm.X.Y.dat](./OUT/SurfaceFluxesNorm.200.1.dat)**

* X : number of depth points
* Y : iteration number
* Columns:

    |  |  |
    | --- | --- |
    1) hnu/kTeff | Dimensionless frequency
    2) Fnu/sigmaTeff^4 | Dimensionless flux
    3) Teff | Effective temperature (constant)
    4) F_logint nu/sigmaTeff^4 | ?

* Source code : [OutputFluxes.c](./OutputFluxes.c)
* [SurfaceFluxesNormt.X.Y.dat](./OUT/SurfaceFluxesNormt.200.1.dat) : same structure, after temperature correction prodcedure

### **[EddingtonFactors.X.Y.dat](./OUT/EddingtonFactors.200.1.dat)**

* X : number of depth points
* Y : iteration number
* Columns:
        
    |  |  |
    | --- | --- |
    1) f| Eddington factor (eq. 20)
    2) h| Eddington factor (eq. 22)
    3) f*J| f * Mean intensity 
    4) h*J| h * Mean intensity 
    5) q| Angle averaged ratio Janist/J 
    6) y| Column depth (g/cm2)
    7) tau| Optical depth
    
* Blocks : One per frequency, X lines each
* Source code : [OutputEddingtonFactors.c](./OutputEddingtonFactors.c)

### **[GauntFactor.dat](./OUT/GauntFactor.dat)**

* Columns:

    |  |  |
    | --- | --- |
    1) logy | log column depth
    2) nu | Frequency (Hz)
    3) T | Temperature (K)
    4) g2 | ?
    5) u | Dimensionless frequency (hnu/kT)
    6) gff | Gaunt factor

* Source code : [McPHAC.c](./McPHAC.c)

### **[InitialkR.dat](./OUT/InitialkR.dat)**

* Columns:

    |  |  |
    | --- | --- |
    1) logy | log column depth
    2) logT | log temperature
    3) rho | density (g/cm3)
    5) P | Pressure (erg/cm3)
    6) kR | ?
    7) kappa | opacity at 1 keV
    8) k | ?

* Source code : [McPHAC.c](McPHAC./.c)

### **[opacities.dat](./OUT/opacities.dat)**

* Columns:

    |  |  |
    | --- | --- |
    1) E | Energy (KeV)
    2) log(kappa) | log opacity at logy = -2
    3) log(kappa) | log opacity at logy = -4
    4) log(kappa) | log opacity at logy = -6

* Blocks : 
* Source code : [PrintOutOpacties.c](./PrintOutOpacties.c)


### **[tau.dat](./OUT/.dat)**

* Columns:

    |  |  |
    | --- | --- |
    1) E | Energy (keV)
    2) y | column depth at tau=1
    3) y | column depth at tau=2
    4) y | column depth at tau=5
    5) y | column depth at tau=10

* Source code : [PrintOutTau.c](./PrintOutTau.c)

### **[Thermo.dat](./OUT/Thermo.dat)**

* Columns:

    |  |  |
    | --- | --- |
    1) logy | Log column depth
    2) logT | Log temperature
    3) P | Pressure (erg/cm3)
    4) rho | density (g/cm3)
    5) kR | ?
 

* Source code : [GetColumnsLog.c](./GetColumnsLog.c)
