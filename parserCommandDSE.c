#include "libDSE.h"
/*This parser works but is more an exercise than something reliable */
void parserCommandDSE(int *argC, char ***argV, char ** filename, double *q0, double *qmax, double *dq, int *ne)
{
	int argc= *argC;
	char **argv= *argV;
	int i, j;
	short int check[4] = {0,0,0,0};
	//printf("%s", argv[1]);
	if((strcmp(argv[1], "-i") == 0)){
			
			*filename= argv[2];
		}
	else {
		puts("no input file!");
		exit(0);
	}
	
	i=2;
	while (i<argc)
	{
		if ((strcmp(argv[i],"-q0") == 0) ){
			*q0=atof(argv[i+1]);
			check[0]=1;
		}
		else {
			if (check[0]==0) *q0=0.0001;
		}
		if ( strcmp(argv[i],"-qmax") ==0){
			*qmax=atof(argv[i+1]);
			check[1]=1;
		}
		else {
			if (check[1]==0) *qmax=15.;
		}
		if (strcmp(argv[i],"-dq") ==0){
			*dq=atof(argv[i+1]);
			check[2]=1;
		}
		else {
			if (check[2]==0) *dq=0.1;
		}
		
		if ((strcmp(argv[i],"-ne")==0)){
			*ne= atoi(argv[i+1]);
			check[3]=1;
		}
		else {
				if (check[3]==0) *ne=1;
		}
		++i;
		//if (i>argc)
			//exit(0);
	}
}
