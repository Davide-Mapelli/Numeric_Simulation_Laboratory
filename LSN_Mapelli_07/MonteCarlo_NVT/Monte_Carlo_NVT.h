/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__

//Random numbers
#include "random.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000; //2
int n_props, iv, iw, igofr;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];


// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres;
// double stima_gofr [m_props], err_gofr[m_props];
double stima_gofr, err_gofr;
double err_pot,err_press;
double deltaV;


//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];

// thermodynamical state
int npart;
double beta,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, norm;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Error2(double,double,int);
void Print(int);

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
