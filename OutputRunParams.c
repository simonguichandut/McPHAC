#include "includes.h"

void OutputRunParams()
// Outputs the input parameters from McPHAC.bash
{
FILE *fp;
fp = fopen("OUT/runparams","w");
fprintf(fp,"Teff        \t%.3e\n"   ,TEFF);
fprintf(fp,"gsurf       \t%.3e\n"   ,GSURFACE);
fprintf(fp,"mincol      \t%.1f\n"   ,MINCOL);
fprintf(fp,"maxcol      \t%.1f\n"   ,MAXCOL);
fprintf(fp,"Ndepths     \t%d\n"     ,NDEPTHS);
fprintf(fp,"Ndepthsnu   \t%d\n"     ,NDEPTHSNU);
fprintf(fp,"Nmu         \t%.d\n"    ,NMU);
fprintf(fp,"Nfreq       \t%.d\n"    ,NFREQ);
fprintf(fp,"Maxfactor   \t%d\n"     ,MAXFACTOR);
fprintf(fp,"Maxfactornu \t%.d\n"    ,MAXFACTORNU);
fprintf(fp,"MaxfracTch  \t%.1e\n"   ,MAXFRACTEMPCHANGE);
fprintf(fp,"Maxcoltau   \t%.3e\n"   ,MAXCOLTAU);
fprintf(fp,"Tguess      \t%.5f\n"   ,TGUESSBDYCOND);
fprintf(fp,"Maxiter     \t%d\n"     ,MAXITER);
fprintf(fp,"anist       \t%d"       ,ANIST);
fclose(fp);
}