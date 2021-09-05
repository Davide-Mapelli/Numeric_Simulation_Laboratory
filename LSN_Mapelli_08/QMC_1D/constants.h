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
This is the definition of the physical constants. The constants have been
put to 1 in order to work with more comfortable numbers (i.e. the potential
energy of the harmonic oscillator on its ground state is 0.25hbar^2*w)
*/

double hbar = 1;
double boltzmann = 1;
double mass = 1;


/****** GLOBAL VARIABLES ************/
/*
PIGS is only a flag variable that is used to determine whether you are running
a zero temperature or a finite temperature Path Integral.
lambda is hbar*hbar/2m, a constant widely used in Path Integral Theory.
dtau is the timestep, that is, the "small imaginary-time" by which the total
propagation time is divided by the Path Integral. Remember that known
propagators are "only" approximations of the true propagator that are valid for
small imaginary-times.
*/

double lambda;
double dtau;
int PIGS;
double alpha;
/*
The following declarations are the variables used by QMC1D. Don't worry, they
are self explaining if you roughly know how a PIMC works.
*/

int timeslices, brownianBridgeReconstructions, brownianBridgeAttempts, brownianMotionReconstructions;
int MCSTEPS, equilibration, blocks, histogram_bins;
int timeslices_averages_start, timeslices_averages_end;
double temperature, imaginaryTimePropagation, delta_variational, delta_translation;
double histogram_start, histogram_end;

int acceptedTranslations, acceptedVariational, acceptedBB, acceptedBM;
int totalTranslations, totalVariational, totalBB, totalBM;

double* positions;
double* potential_energy;
double* potential_energy_accumulator;
double* potential_energy_square_accumulator;
                                                                                                                 
double* kinetic_energy;
double* kinetic_energy_accumulator;
double* kinetic_energy_square_accumulator;
                                                                                                                 
double* positions_histogram;
double* positions_histogram_accumulator;
double* positions_histogram_square_accumulator;

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
