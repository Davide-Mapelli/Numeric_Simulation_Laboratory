/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/*
The comments in this file are not detailed. See qmc1d.cpp for a better explanation.
*/

void readInput();   // reads input from the file "input.dat"
void deleteMemory(); // handles the dynamic allocation of memory
void initialize();  // initializes the variables
void consoleOutput(); // writes the output on the screen
                                                                                                                 
                                                                                                                 
double potential_density_matrix(double val, double val_next);
double u_prime(double x, int m);
double u_sec(double x, int m);

/*
potential_density_matrix returns only the potential part of the correlation between two adjacent timeslices.
*/

double external_potential(double);  // this is the external potential definition
double external_potential_prime(double); // ...and here goes its first derivative
double external_potential_second(double); // ... and its second derivative 

/*
The derivatives are necessary for the evaluation of the kinetic estimator, because it contains
the laplacian operator ! 
*/                                                                                                    
                                                                                                                 
void translation(); // performs a rigid translation
void brownianBridge();  // reconstructs a segment of the polymer with a free particle propagation. 
void brownianMotion(int);  // reconstructs a segment at the extremities of the polymer with a free particle propagation. 
                                                                                                                 
double variationalWaveFunction(double);  
/*variationalWaveFunction is the variational wave function that is
projected in a PIGS simulation.
*/
double variationalWaveFunction_second(double);
double variationalLocalEnergy(double val);
/*
as for the potential, you have to specify its first and second derivative for the evaluation
of the kinetic local energy.
*/                                                                                                                

int index_mask(int); 
/* index_mask is just a compatibility function that takes into account whether the
polymer is open (PIGS) or closed in periodic boundary contitions (PIMC-ring polymer).
*/

void upgradeAverages(); // at every MCSTEP accumulates the estimators values.

void upgradeHistogram(); // fills the histogram of positions foreach MCSTEP
void endBlock(); // finalizes the averages at the end of each block

double kineticEstimator(double,double);  // evaluates the kinetic energy along the polymer
void finalizePotentialEstimator();
void finalizeKineticEstimator();
void finalizeHistogram();
/*
The last three functions are called at the end of the simulation, basically they average over each
block and evaluate the error on the block average. This is an application of the central limit
theorem.
*/

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
