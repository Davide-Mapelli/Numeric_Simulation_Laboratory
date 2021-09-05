/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/*********************** QMC1D **************************
************** PATH INTEGRAL GROUND STATE ***************
************** PATH INTEGRAL MONTE CARLO ****************
************ APPLIED TO A SINGLE PARTICLE ***************
************** IN AN EXTERNAL POTENTIAL ****************/
/*
NOTE: you need the root package to be installed before compiling this program.
See: http://root.cern.ch
There are two other source files, too:
constants.h: contains the declaration of every global variable that has been used.
functions.h: contains the declaration of the function with a brief description.
Once compiled, QMC1D is invoked with the command: "./qmc1d". It will read the settings 
in the file "input.dat".
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <TRandom3.h>
#include "constants.h"
#include "functions.h"

#define LEFT 0
#define RIGHT 1

TRandom3* generator;

using namespace std;

int main()
{
        readInput();  
	initialize();  
/* at this time, every variable you see, such for instance "equilibration",
has been either acquired from "input.dat" by the readInput() function or
opportunely initialized by the initialize() function. */
	for(int i=0;i<equilibration;i++)
	{
		if(PIGS)   // only a PIGS polymer has a start and an end. 
		{
			brownianMotion(LEFT);
			brownianMotion(RIGHT);
		}
		
		translation();
		
		for(int j=0;j<brownianBridgeAttempts;j++) 
			brownianBridge();
	}
	
	for(int b=0;b<blocks;b++)
	{
		for(int i=0;i<MCSTEPS;i++)
		{
			if(PIGS)
			{
				brownianMotion(LEFT);
				brownianMotion(RIGHT);
			}
			translation();
		
			for(int j=0;j<brownianBridgeAttempts;j++)
				brownianBridge();
			
			upgradeAverages();
		}
		cout<<"Completed block: "<<b+1<<"/"<<blocks<<endl;
		endBlock();
	}
	
	consoleOutput();
	finalizePotentialEstimator();
	finalizeKineticEstimator();
	finalizeHistogram();

	deleteMemory();  // de-allocate dynamic variables.
	return 0;
}

// This is the primitive approximation without the kinetic correlation
double potential_density_matrix(double val, double val_next)
{
	double dens_left = -dtau*external_potential(val)/2;
	double dens_right = -dtau*external_potential(val_next)/2;
	
	return dens_left+dens_right;
}

// Initialization of the variables and allocation of the memory.
void initialize()
{
	lambda = hbar*hbar/(2*mass);
	if(temperature==0)
		PIGS=1;
	else
		PIGS=0;
	
	if(PIGS)
		dtau = imaginaryTimePropagation/(timeslices-1);
	else
		dtau = hbar/(boltzmann*temperature*timeslices);
	
	acceptedTranslations=0;
	acceptedVariational=0;
	acceptedBB=0;
        acceptedBM=0;
	totalTranslations=0;
	totalVariational=0;
	totalBB=0;
        totalBM=0;
	
	generator = new TRandom3();
	
	positions=new double[timeslices];
	potential_energy=new double[timeslices];
	potential_energy_accumulator=new double[timeslices];
	potential_energy_square_accumulator=new double[timeslices];
                                                                                                                 
	kinetic_energy=new double[timeslices];
	kinetic_energy_accumulator=new double[timeslices];
	kinetic_energy_square_accumulator=new double[timeslices];
                                                                                                                
	positions_histogram=new double[histogram_bins];
	positions_histogram_accumulator=new double[histogram_bins];
	positions_histogram_square_accumulator=new double[histogram_bins];
	
	for(int i=0;i<timeslices;i++)
		positions[i]=0.0;
	
	for(int i=0;i<timeslices;i++)
	{
		potential_energy[i]=0;
		potential_energy_accumulator[i]=0;
		potential_energy_square_accumulator[i]=0;

		kinetic_energy[i]=0;
		kinetic_energy_accumulator[i]=0;
		kinetic_energy_square_accumulator[i]=0;
	}
	
	for(int i=0;i<histogram_bins;i++)
	{
		positions_histogram[i]=0;
		positions_histogram_accumulator[i]=0;
		positions_histogram_square_accumulator[i]=0;
	}
	alpha=0;
}

// The external potential. You can modify this function but don't forget
// to modify its first and second derivatives too !
double external_potential(double val)
{
	double k_elastic = 1;
	return k_elastic*val*val/2.0;
}

double external_potential_prime(double val)
{
	double k_elastic =1;
	return k_elastic*val;
}

double external_potential_second(double val)
{
	double k_elastic=1;
	return k_elastic;
}

// The same applies to the variational Wave Function...
// You can modify this function but don't forget
// to modify its second derivative below!
double variationalWaveFunction(double v)
{
	return 1.0;
	//return exp(-0.5*v*v);
}

double variationalWaveFunction_second(double v)
{
	return 0;
	//return v*v*exp(-0.5*v*v) - exp(-0.5*v*v);
}

void translation()
{
	totalTranslations++;
	double delta = generator->Uniform(-delta_translation,delta_translation);
	double acc_density_matrix_difference=0;
	int last = timeslices;
	if(PIGS)
		last=timeslices-1;
		
	for(int i=0;i<last;i++)
	{
		int inext = index_mask(i+1);
		double newcorr,oldcorr;
		newcorr = potential_density_matrix(positions[i]+delta,positions[inext]+delta);
		oldcorr = potential_density_matrix(positions[i],positions[inext]);
		acc_density_matrix_difference += oldcorr-newcorr;
	}
	// metropolis: PIGS contains also the statistical weight of the variational Wave Function.
	double acceptance_probability = exp(-acc_density_matrix_difference);
	
	if(PIGS)
		acceptance_probability *=variationalWaveFunction(positions[0]+delta)*variationalWaveFunction(positions[timeslices-1]+delta)/(variationalWaveFunction(positions[0])*variationalWaveFunction(positions[timeslices-1]));
	
	if(generator->Rndm()<acceptance_probability)
	{
		for(int i=0;i<timeslices;i++)
			positions[i]+=delta;
		acceptedTranslations++;
	}
}

/* BB removes a segment of the polymer, in this case from "starting_point+1" to "endpoint-1"
and replaces it with a free particle propagation. The free particle propagation is achieved
with the gaussian sampling of the kinetic part of the density matrix.
The function index_mask handles the compatibility between PIGS and PIMC: in PIGS the polymer is 
open, so you can't have a starting index greater than an ending index. In PIMC, instead, you
have a ring polymer so when you reach the end you can continue from the beginning. 
The compatibility solution that has been chosen consists in viewing the ring polymer as an open
polymer that has been closed on periodic boundary contitions. index_mask takes this into account. */
void brownianBridge()
{
	totalBB++;
	int available_starting_points = timeslices-brownianBridgeReconstructions-1; // for PIGS simulation
	if(!PIGS)
		available_starting_points = timeslices-1;
	int starting_point = (int)(generator->Rndm()*available_starting_points);
	
	int endpoint = index_mask(starting_point + brownianBridgeReconstructions + 1);
	
	double starting_coord = positions[starting_point];
	double ending_coord = positions[endpoint];
	double new_segment[brownianBridgeReconstructions+2];
	new_segment[0]=starting_coord;
	new_segment[brownianBridgeReconstructions+1]=ending_coord;
	double previous_position = starting_coord;
	for(int i=0;i<brownianBridgeReconstructions;i++)
	{
		int left_reco = brownianBridgeReconstructions-i;
		// gaussian sampling of the free particle propagator
		double average_position = previous_position + (ending_coord-previous_position)/(left_reco+1);
		double variance = 2*lambda*dtau*left_reco/(left_reco+1);
		double newcoordinate = generator->Gaus(average_position,sqrt(variance));
		new_segment[i+1] = newcoordinate;
		previous_position=newcoordinate;
	}
	
	// metropolis. Note that the kinetic part has been sampled exactely, thus only the
	// potential part of the density matrix determines the acceptance probability of the move.
	double acc_density_matrix_difference=0;
	for(int i=0;i<brownianBridgeReconstructions+1;i++)
	{
		int i_old = index_mask(starting_point+i);
		int i_next_old = index_mask(starting_point+i+1);
		double newcorr,oldcorr;
		newcorr = potential_density_matrix(new_segment[i],new_segment[i+1]);
		oldcorr = potential_density_matrix(positions[i_old],positions[i_next_old]);
		acc_density_matrix_difference += oldcorr-newcorr;
	}
	
	double acceptance_probability = exp(-acc_density_matrix_difference);
	if(generator->Rndm()<acceptance_probability)
	{
		for(int i=1;i<brownianBridgeReconstructions+1;i++)
		{
			int i_old = index_mask(starting_point+i);
			positions[i_old]=new_segment[i];
		}
		acceptedBB++;
	}
}

/* BM removes a segment at one of the two ends of the polymer,
and replaces it with a free particle propagation using a Brownian Bridge after the sampling of the
starting (left move) or final (right move) position. The free particle propagation is achieved
with the gaussian sampling of the kinetic part of the density matrix. */
void brownianMotion(int which) // BM is called only for PIGS simulations
{
	int starting_point, endpoint, left_reco;
        double starting_coord, ending_coord, average_position, variance, newposition, oldposition;

        totalBM++;

        if(which==LEFT)
        {
                starting_point = 0;
                endpoint = brownianMotionReconstructions+1;
		ending_coord = positions[endpoint];
		average_position = ending_coord;
                variance = 2*lambda*dtau*(brownianMotionReconstructions+1);
		starting_coord = generator->Gaus(average_position,sqrt(variance));
                oldposition = positions[starting_point];
                newposition = starting_coord;
        }
        else
        {
		starting_point = timeslices-2-brownianMotionReconstructions;
		endpoint = timeslices-1;
		starting_coord = positions[starting_point];
                average_position = starting_coord;
                variance = 2*lambda*dtau*(brownianMotionReconstructions+1);
		ending_coord = generator->Gaus(average_position,sqrt(variance));
		oldposition = positions[endpoint];
		newposition = ending_coord;
        }

        double new_segment[brownianMotionReconstructions+2];
        new_segment[0]=starting_coord;
        new_segment[brownianMotionReconstructions+1]=ending_coord;
        double previous_position = starting_coord;
        for(int i=0; i<brownianMotionReconstructions; i++)
        {
                left_reco = brownianMotionReconstructions-i;
                // gaussian sampling of the free particle propagator
                average_position = previous_position + (ending_coord-previous_position)/(left_reco+1);
                variance = 2*lambda*dtau*left_reco/(left_reco+1);
                double newcoordinate = generator->Gaus(average_position,sqrt(variance));
                new_segment[i+1] = newcoordinate;
                previous_position=newcoordinate;
        }

        // metropolis. Note that the kinetic part has been sampled exactely, thus only the
        // potential part of the density matrix determines the acceptance probability of the move.
        double acc_density_matrix_difference=0;
        for(int i=0;i<brownianMotionReconstructions+1;i++)
        {
                double newcorr,oldcorr;
                newcorr = potential_density_matrix(new_segment[i],new_segment[i+1]);
                oldcorr = potential_density_matrix(positions[starting_point+i],positions[starting_point+i+1]);
                acc_density_matrix_difference += oldcorr-newcorr;
        }

        double acceptance_probability = exp(-acc_density_matrix_difference)*variationalWaveFunction(newposition)/variationalWaveFunction(oldposition);
        if(generator->Rndm()<acceptance_probability)
        {
                for(int i=0;i<brownianMotionReconstructions+2;i++)
                {
                        positions[starting_point+i]=new_segment[i];
                }
                acceptedBM++;
        }
}

int index_mask(int ind)
{
	if(PIGS)
		return ind;  // no pbc over indices
	else
	{
		int new_ind=ind;
		while(new_ind>=timeslices)   // pbc over indices
			new_ind-=timeslices;
		return new_ind;
	}
}

void consoleOutput()
{
	cout<<"Acceptances:"<<endl;
	if(PIGS)
		cout<<"BM: "<<((double)acceptedBM)/totalBM<<endl;
	cout<<"Transl: "<<((double)acceptedTranslations)/totalTranslations<<endl;
	cout<<"BB: "<<((double)acceptedBB)/totalBB<<endl;
}


/* This function accumulates the expectation values in their respective variables. 
   At the end of the block, these variables are divided by the MCSTEPS value and the block average
   and its squared value are accumulated in apposite variables.  
 */
void upgradeAverages()
{
	for(int i=0;i<timeslices;i++)
		potential_energy[i]+=external_potential(positions[i]);
	
	int flag=0;
	if(PIGS)
		flag=1;
		
	for(int i=flag;i<timeslices-flag;i++)
	{
		int i_mod = index_mask(i);
		int i_mod_next = index_mask(i+1);
		kinetic_energy[i_mod]+= kineticEstimator(positions[i_mod],positions[i_mod_next]);
	}
	if(flag)
	{
		kinetic_energy[0]+=variationalLocalEnergy(positions[0]);
		kinetic_energy[timeslices-1]+=variationalLocalEnergy(positions[timeslices-1]);
	}
	
	upgradeHistogram();
}

/*
This functions performs a common way to fill in the values of an histogram. 
It's quite straightforward.
*/
void upgradeHistogram()
{
	double delta_pos = (histogram_end-histogram_start)/histogram_bins;
	for(int i=timeslices_averages_start; i<=timeslices_averages_end; i+=1)
	{
		int k=1;
		while((histogram_start + k*delta_pos)<positions[i])
			k++;
		positions_histogram[k-1]+=1;
	}
}

void endBlock()  // calculating and accumulating block averages
{
	for(int i=0;i<timeslices;i++)
	{
		potential_energy[i]/=MCSTEPS;
		potential_energy_accumulator[i]+=potential_energy[i];
		potential_energy_square_accumulator[i]+=potential_energy[i]*potential_energy[i];
		potential_energy[i]=0;
		kinetic_energy[i]/=MCSTEPS;
		kinetic_energy_accumulator[i]+=kinetic_energy[i];
		kinetic_energy_square_accumulator[i]+=kinetic_energy[i]*kinetic_energy[i];
		kinetic_energy[i]=0;
		
	}
	
	for(int i=0;i<histogram_bins;i++)
	{
		positions_histogram[i]/=MCSTEPS;
		positions_histogram_accumulator[i]+=positions_histogram[i];
		positions_histogram_square_accumulator[i]+=positions_histogram[i]*positions_histogram[i];
		positions_histogram[i]=0;
	}
}

/*
finalize**** functions average the accumulators of the block averages and their square values.
Then the error is calculated by the usual formula err(A)=sqrt(|<A><A>-<A*A>|/Nblocks)
*/
void finalizePotentialEstimator()
{
	ofstream out("potential.dat");
	for(int i=0;i<timeslices;i++)
	{
		double potential_energy_average = potential_energy_accumulator[i]/blocks;
		double potential_energy_square_avg = potential_energy_square_accumulator[i]/blocks;
		double p_error =sqrt(abs(potential_energy_average*potential_energy_average-potential_energy_square_avg)/blocks);
		out<<i<<" "<<potential_energy_average<<" "<<p_error<<endl;
	}
	out.close();
}

void finalizeKineticEstimator()
{
	ofstream out("kinetic.dat");
	for(int i=0;i<timeslices;i++)
	{
		double kinetic_energy_average = kinetic_energy_accumulator[i]/blocks;
		double kinetic_energy_square_avg = kinetic_energy_square_accumulator[i]/blocks;
		double k_error =sqrt(abs(kinetic_energy_average*kinetic_energy_average-kinetic_energy_square_avg)/blocks);
		out<<i<<" "<<kinetic_energy_average<<" "<<k_error<<endl;

	}
	out.close();
}

void finalizeHistogram()
{
	ofstream out("probability.dat");
        double current_position, hist_average, hist_square_avg, error;
	double delta_pos = (histogram_end-histogram_start)/histogram_bins;
        double norma = 0.0;
	for(int i=0; i<histogram_bins; i++)
	{
		norma += positions_histogram_accumulator[i]/blocks;
	}
        norma *= delta_pos;
        for(int i=0; i<histogram_bins; i++)
        {
                current_position = histogram_start + (i+0.5)*delta_pos;
		hist_average = positions_histogram_accumulator[i]/blocks;
		hist_square_avg = positions_histogram_square_accumulator[i]/blocks;
		error =sqrt(abs(hist_average*hist_average-hist_square_avg)/blocks);
                out << current_position << " " << hist_average/norma << " " << error/norma << endl;
                
        }
	out.close();
}

// (-hbar*hbar/2m)d^2/dx^2G(x,x',dtau)
double kineticEstimator(double value,double next_value)
{
	double kinetic_prime = (value-next_value)/(2*lambda*dtau);
	double kinetic_second= 1./(2*lambda*dtau);
	double term_1 = (dtau/2)*external_potential_prime(value)+kinetic_prime;
	double term_2 = (dtau/2)*external_potential_second(value)+kinetic_second;
	return -(hbar*hbar/(2*mass))*(term_1*term_1 - term_2);
}

// (-hbar*hbar/2m)(d^2/dx^2G(x,x',dtau))/G(x,x',dtau)
double variationalLocalEnergy(double val)
{
	double psi = variationalWaveFunction(val);
	double laplacian_psi = variationalWaveFunction_second(val);
	return -(hbar*hbar/(2*mass))*laplacian_psi/psi;
}

void readInput()
{

	ifstream input_file("input.dat");
	char* string_away = new char[60];

	input_file >> string_away >> timeslices;
	input_file >> string_away >> temperature;
	input_file >> string_away >> imaginaryTimePropagation;
	input_file >> string_away >> brownianMotionReconstructions;
	input_file >> string_away >> delta_translation;
	input_file >> string_away >> brownianBridgeReconstructions;
	input_file >> string_away >> brownianBridgeAttempts;
	input_file >> string_away >> MCSTEPS;
	input_file >> string_away >> equilibration;
	input_file >> string_away >> blocks;
	input_file >> string_away >> histogram_bins;
	input_file >> string_away >> histogram_start;
	input_file >> string_away >> histogram_end;
	input_file >> string_away >> timeslices_averages_start>>timeslices_averages_end;
	input_file.close();
	delete [] string_away;
}

void deleteMemory()
{
	delete positions;
	delete potential_energy;
	delete potential_energy_accumulator;
	delete potential_energy_square_accumulator;
                                                                                                                 
	delete kinetic_energy;
	delete kinetic_energy_accumulator;
	delete kinetic_energy_square_accumulator;
                                                                                                                 
	delete positions_histogram;
	delete positions_histogram_accumulator;
	delete positions_histogram_square_accumulator;

        delete generator;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
