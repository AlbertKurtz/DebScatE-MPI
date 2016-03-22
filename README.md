# DebScatE-MPI
Code that calculates the Debye Scattering Equation for a given set of atoms (4 columns file) using MPI.

Calculates the scattered intensity given a range for the scattering vector q= 4*Pi/(lambda) *sin(theta).
I(q) = Sum_{i=0}^{N}Sum_{j=0}^{N} f_j(q)f_i(q)Sinc(q*r_{ij})

Once compiled it works like this: mpirun -np n ./a.out -i "inputfile" -q0 "q0" -qmax "qmax" -dq "dq" -ne "ne"
where
n is number of threads; 
q0 the minimum value for the scattering vector q;
qmax the maximum value for the scattering vector q;
dq is the step of q;
ne is the number of atomic elements in the sample (NOTE: code works only for one species, needs work).
