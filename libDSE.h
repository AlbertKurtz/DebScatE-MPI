#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <mpi.h>
#include <string.h>

#define FFNUMPAR 9
#define mpi_root 0

#define dist(x1, y1, z1, x2, y2, z2) sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
//#define dw(q,dwpar) exp(-dwpar*q*q/(8*M_PI*M_PI))
//#define ff(q, ffpar) (ffpar[0]*exp(-ffpar[1]*q*q/(16*M_PI*M_PI)))+(ffpar[2]*exp(-ffpar[3]*q*q/(16*M_PI*M_PI)))+(ffpar[4]*exp(-ffpar[5]*q*q/(16*M_PI*M_PI)))+(ffpar[6]*exp(-ffpar[7]*q*q/(16*M_PI*M_PI)))+ffpar[8]



void  parserxyzDSE(char *fileneeded, double **atomx, double **atomy, double **atomz, char (**atomtype)[5], unsigned int *N);
void parserCommandDSE(int *argC, char ***argV, char ** filename, double *q0, double *qmax, double *dq, int *ne);
void ffparamDSE(char *atype, double (*ffpar)[FFNUMPAR], double *dwpar);
double ff(double q, double *ffpar);


