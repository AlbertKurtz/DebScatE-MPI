#include "libDSE.h"

void ffparamDSE(char *atype, double (*ffpar)[FFNUMPAR], double *dwpar)
{
	FILE *fp, *gp;
	fp = fopen ("formfactorsDSE", "r");
	gp = fopen ("dwparamDSE", "r");
	
	int i;
	char type[strlen(atype)], type2[strlen(atype)];
	
	for (i=0; i<(strlen(atype)+1); ++i)
	{
		atype[i]= toupper(atype[i]);
	}
	
	for (i=0; i< 300; ++i)
	{
		fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf", type, &(*ffpar)[0], &(*ffpar)[1], &(*ffpar)[2], &(*ffpar)[3], &(*ffpar)[4], &(*ffpar)[5], &(*ffpar)[6], &(*ffpar)[7], &(*ffpar)[8] );
		fscanf(gp,"%s %lf", type2, dwpar);
		if (strcmp(atype,type) == 0)
		{
			i=301;
		}
	}
	
	fclose (fp);
	fclose (gp);
	
}

double ff(double q, double *ffpar)
{
	return (ffpar[0]*exp(-ffpar[1]*q*q/(16*M_PI*M_PI)))+(ffpar[2]*exp(-ffpar[3]*q*q/(16*M_PI*M_PI)))+(ffpar[4]*exp(-ffpar[5]*q*q/(16*M_PI*M_PI)))+(ffpar[6]*exp(-ffpar[7]*q*q/(16*M_PI*M_PI)))+ffpar[8];
}

#define eightpipi 8*M_PI*M_PI

double dw(double q, double dwpar)
{
	return exp(-dwpar*q*q/eightpipi);
}



