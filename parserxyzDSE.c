#include "libDSE.h"

//NOTE:  the sensitive point is the string length, in this case is 5 but can change


void  parserxyzDSE(char *fileneeded, double **atomx, double **atomy, double **atomz, char (**atomtype)[5], unsigned int *N)
{
	FILE * fp;	
	unsigned int i;
	
	fp = fopen(fileneeded, "r");
	if (fp == NULL)	{
	  exit(EXIT_FAILURE);
	}
	fscanf(fp,"%u", N);
	
	*atomx = (double*) malloc(sizeof(double)*(*N));
	*atomy = (double*) malloc(sizeof(double)*(*N));
	*atomz = (double*) malloc(sizeof(double)*(*N));
	*atomtype= malloc(5*sizeof(char) * (*N));

	for (i=0; i<(*N); ++i)
	{		
		
		fscanf(fp,"%s %lf %lf %lf", (*atomtype)[i],  &(*atomx)[i], &(*atomy)[i], &(*atomz)[i] ); 
	}	
	

	fclose(fp);
}
