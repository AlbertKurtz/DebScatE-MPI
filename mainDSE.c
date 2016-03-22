#include "libDSE.h"

int numnodes, myid;

int main (int argc, char* argv[])
{

	MPI_Init(&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &numnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);	
	
	unsigned int N, i, j, k, Nsteps, m;
	double *ax, *ay, *az, *Iplotpar, *Iplot;
	double q0, qmax, dq, rij, Q;
	int ne;
	char (*at)[5];
	char * datafile;
	double ffpar[FFNUMPAR];
	double dwpar;
	double time0, time1;
	
	time0=MPI_Wtime();
	
	if (myid==mpi_root)
	{
		
		parserCommandDSE( &argc, &argv, &datafile, &q0, &qmax, &dq, &ne);
		parserxyzDSE(datafile, &ax, &ay, &az, &at, &N);	
		ffparamDSE(at[0], &ffpar, &dwpar);
		Nsteps= (qmax-q0)/dq;
		
	}
	
	MPI_Bcast (&N, 1, MPI_UNSIGNED, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (&q0, 1, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (&qmax, 1, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (&dq, 1, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (&Nsteps, 1, MPI_UNSIGNED, mpi_root, MPI_COMM_WORLD);
	
	//Qplot = (double*) malloc(sizeof(double)*(Nsteps));
	Iplotpar = (double*) malloc(sizeof(double)*(Nsteps));
	Iplot = (double*) malloc(sizeof(double)*(Nsteps));
	
	m= (unsigned int) N/numnodes;
	
	if(myid!=mpi_root)
	{
		ax = (double*) malloc(sizeof(double)*(N));
		ay = (double*) malloc(sizeof(double)*(N));
		az = (double*) malloc(sizeof(double)*(N));
		at= malloc(5*sizeof(char) * (N));
	}
	
	MPI_Bcast (ax, N, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (ay, N, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (az, N, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (at, 5*N, MPI_CHAR, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (ffpar, FFNUMPAR, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	MPI_Bcast (&dwpar, 1, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
	
	
	//one of the sums the sum is split into m which is N/numnodes
	
	if (myid != numnodes-1)
	{
		for(k=0; k<Nsteps; ++k)
		{
			Iplotpar[k]=0.;
			Q=k*dq;
			for (i=myid*m; i<(myid+1)*m; ++i)
			{
				for (j=0; j<N; ++j)
				{
					rij=dist(ax[i], ay[i], az[i], ax[j], ay[j], az[j]);
					rij*=Q;
					rij+=0.0000001;
					Iplotpar[k]+=sin(rij)/(rij);
					//printf("%lf\n", Iplotpar[k]);
				}
			}
			
			Iplotpar[k]*= ff(Q, ffpar)* ff(Q, ffpar); //*dw(Q, dwpar);
			//printf("%lf\n", Iplotpar[k]);
		
			
		}
	}
	
	if (myid== numnodes-1)
	{
		for(k=0; k<Nsteps; ++k)
		{
			Iplotpar[k]=0.;
			Q=k*dq;
			for (i=myid*m; i<N; ++i)
			{
				for (j=0; j<N; ++j)
				{
					rij=dist(ax[i], ay[i], az[i], ax[j], ay[j], az[j]);
					rij*=Q;
					rij+=0.0000001;
					Iplotpar[k]+=sin(rij)/(rij);
					//printf("%lf\n", Iplotpar[k]);
				}
			}
			
			Iplotpar[k]*= ff(Q, ffpar)* ff(Q, ffpar); //*dw(Q, dwpar);
			//printf("%lf\n", Iplotpar[k]);
		
			
		}
	}
	free(at);
	free(ax);
	free(ay);
	free(az);
	
	MPI_Reduce (Iplotpar, Iplot, Nsteps, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);	
	free(Iplotpar);
	if (myid==mpi_root)
	{
		/*
		double  Qplot[Nsteps];
		for (k=0; k<Nsteps; ++k)
		{
			Qplot[k]=dq*k;
		}
		*/
		
		FILE *fp;
		fp= fopen("intensity.tsv", "w+");
		for (k=0; k<Nsteps; ++k)
		{
			fprintf(fp,"%lf\t%lf\n", dq*k, Iplot[k]);
			printf("%lf\t%lf\n", dq*k, Iplot[k]);	
		}
		
		fclose(fp);
	}
	time1=MPI_Wtime();
	if(myid==mpi_root)
		printf("\ntime elapsed %lf\n", time1-time0);
	
	MPI_Finalize();
}

