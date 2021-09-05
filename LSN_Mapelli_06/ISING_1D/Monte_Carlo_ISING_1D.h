/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __ISING__
#define __ISING__

//Random numbers
#include "random.h"
#include <vector>
#include <cstdlib>

using namespace std;
int seed[4];
Random rnd;

//parameters, observables
//const int m_props=1000;
const int n_props=4;
int iu,ic,im,ix,ig;
double nbins;
double walker[n_props]; //m?

// averages
double blk_av[n_props], blk_norm, accepted, attempted; //m
double glob_av[n_props],glob_av2[n_props];
double stima_u,stima_c,stima_m,stima_x,stima_g;
double err_u,err_c,err_m,err_x,err_g;


//configuration
const int m_spin=50;
double s[m_spin];

// thermodynamical state
int nspin;
double beta,temp,J,h;
double range = 2-0.5;
int step = 10;


// simulation
int nstep, nblk, metro;
bool type;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);
bool search(int, vector <int>);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
