/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props, igofr;
double bin_size;
const int nbins=100;
int bin[100];
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stimagofr;
double deltaV;
// averages
double acc,att;
double blk_av[nbins],blk_norm,accepted,attempted;
double glob_av[nbins],glob_av2[nbins];
double stima_gofr [nbins], err_gofr[nbins];
double norm;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,tempstar,vol,rho,box,rcut;

// simulation
bool type;
int nstep, iprint, seed;
double delta;

//functions
int funzione(void);
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfFinal(bool type);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Average(void);
void Gofr(void);
double error (double, double, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
