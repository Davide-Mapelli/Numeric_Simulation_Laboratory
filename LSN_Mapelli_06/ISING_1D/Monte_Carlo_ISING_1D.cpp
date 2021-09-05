/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include "Monte_Carlo_ISING_1D.h"
//#include "random.h"

using namespace std;

int main(){ 
  h = 0;
  for (int w = 0; w < 2; w++){
  metro = 0;
  for (int m = 0; m<2; m++){
    temp = 0.5;
    for (int t = 0; t<10; t++){
      Input(); //Inizialization
      cout << "initial configuration : " << endl;
      for (int n = 0; n < nspin; n++)cout << s[n] << "  ";
      cout << endl;
      for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep){
          Move(metro);
          Measure();
          Accumulate(); //Update block averages
        } 
        Averages(iblk);   //Print results for current block
      }
      ConfFinal(); //Write final configuration
      temp += range / step;
    }
    metro += 1;
  }
  h += 0.02;
  }
  return 0;
}
///////////////////////////////////////////////////////////////////
////////////////////////FUNCTIONS//////////////////////////////////
///////////////////////////////////////////////////////////////////

////////////////////////INPUT//////////////////////////////////

void Input(void){
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");
  ReadInput >> type;
  //ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;
  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;
  ReadInput >> J;
  cout << "Exchange interaction (J)= " << J << endl;
  //ReadInput >> h;
  cout << "External field (h) = " << h << endl << endl;
  //ReadInput >> metro; // if(metro=0) {Metropolis}; else {Gibbs}
  ReadInput >> nblk;
  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  //n_props = 4; //Number of observables

//initial configuration
  if (type == 0 ){
    cout << "Starting from a casual configuration" << endl << endl;
   for (int i=0; i<nspin; ++i){
      if(rnd.Rannyu() >= 0.5) s[i] = +1;
      else s[i] = -1;
      //cout << "s[" << i << "] = " << s[i] << endl; 
    }
  }
  else{
    cout << "Starting from the previous configuration" << endl << endl;
    ReadInput.open("config.final");
    for (int i=0; i<nspin; ++i) ReadInput >> s[i] ; 
    ReadInput.close();
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();
//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  //cout << "Initial heat capacity = " << walker[ic]/(double)nspin << endl;
  //cout << "Initial magnetization = " << walker[im]/(double)nspin << endl;
  //cout << "Initial magnetic susceptibility = " << walker[ix]/(double)nspin << endl;
}

///////////////////////////////////////////////////////////

////////////////////////MOVE//////////////////////////////////

void Move(int metro){
  int o;
  double p, energy_old, energy_new, sm, deltaE;
  //double energy_up, energy_down;
  vector <int> saved;
  
  for(int i=0; i<nspin; i++){ //++i
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    do{o = (int)(rnd.Rannyu()*nspin);}
    while (search( o, saved ));
    saved.push_back(o);
    sm = s[o];
    energy_old = Boltzmann(sm,o);
    energy_new = Boltzmann(-sm,o);
    deltaE = energy_new - energy_old;

    if(metro==0){ //Metropolis
      //cout << "   deltaE =  " << deltaE;
     
      if ( deltaE < 0){
        s[o] *= -1;
        accepted++;
        //cout << accepted << endl;
        }
      else{
        p = rnd.Rannyu();
        if (p < exp(-beta*deltaE)){
          s[o] *= -1;
          accepted ++;
          //cout << accepted << endl;
        }
      }
    }
    else{ //Gibbs sampling
    double x = rnd.Rannyu();
    double p = 1.0 / (1.0 + exp( beta*deltaE) );
      if (x < p) {
        s[o] *= -1;
        accepted++;
      }
    }
    
  saved.clear();
  }
}

//SEARCH
bool search (int k, vector <int> saved){
  for (unsigned int i =0; i< saved.size(); i++)if (saved[i] == k)return true;
  return false;
}
/////////////////////////////


////////////////////////BOLTZMANN//////////////////////////////////

double Boltzmann(int sm, int ip){
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

////////////////////////MEASURE//////////////////////////////////

void Measure(){
  //int bin;
  double u = 0.0, m = 0.0; 
  //double u2 = 0.0, m2 = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]); //<H>
    //u2 += u*u;  //<H^2>
    m += s[i];  //<M>
    //m2 *= m*m;  //<M^2>
  }
// INCLUDE YOUR CODE HERE
  
  walker[iu] = u ;
  walker[ic] = u*u ;
  walker[im] = m;
  walker[ix] = m*m;
// INCLUDE YOUR CODE HERE
}

////////////////////////RESET//////////////////////////////////

void Reset(int iblk) //Reset block averages
{
  if(iblk == 1){
    for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i) blk_av[i] = 0; 
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

////////////////////////ACCUMULATE//////////////////////////////////

void Accumulate(void) //Update block averages
{
  for(int i=0; i<n_props; ++i)blk_av[i] += walker[i];
  blk_norm += 1.0;
}

////////////////////////AVERAGES//////////////////////////////////


void Averages(int iblk) //Print results for current block
{
  ofstream Ene, Heat, Mag, Chi;
  ofstream AveE, AveH, AveM, AveX;
  //const int wd=12;
  attempted = nspin * nstep;
    
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
  
  /////Energy
  Ene.open("output.ene."+to_string(metro),ios::app);
  stima_u = blk_av[iu]/blk_norm/(double)nspin; 
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);
  //Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
  if ((temp == 0.5) & (h == 0)) Ene << glob_av[iu]/(double)iblk << " " << err_u << " ";
  AveE.open("ave.ene."+to_string(metro),ios::app);
  if ((iblk == nblk) & (h == 0))  AveE << glob_av[iu]/(double)iblk << " " << err_u << " ";
  AveE.close();
  Ene.close();

  /////Heat capacity
  Heat.open("output.heat."+to_string(metro),ios::app);
  stima_c = beta*beta*(blk_av[ic]/  blk_norm / (double)nspin - stima_u*stima_u*(double)nspin); 
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);
  //Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  if ((temp == 0.5) & (h == 0)) Heat << glob_av[ic]/(double)iblk << " " << err_c << " ";
  AveH.open("ave.heat."+to_string(metro),ios::app);
  if ((iblk == nblk) & (h == 0))  AveH << glob_av[ic]/(double)iblk << " " << err_c << " ";
  AveH.close();
  Heat.close();

  /////Magnetization
  Mag.open("output.mag."+to_string(metro),ios::app);
  stima_m = blk_av[im]/blk_norm/(double)nspin; 
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);
  //Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
  if ((temp == 0.5) & (h == 0.02)) Mag << glob_av[im]/(double)iblk << " " << err_m << " ";
  AveM.open("ave.mag."+to_string(metro),ios::app);
  if ((iblk == nblk) & (h == 0.02))  AveM << glob_av[im]/(double)iblk << " " << err_m << " ";
  AveM.close();
  Mag.close();

  /////Magnetic Susceptibility (X)
  Chi.open("output.chi."+to_string(metro),ios::app);
  stima_x = beta * blk_av[ix]/blk_norm/(double)nspin; 
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);
  //Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  if ((temp == 0.5) & (h == 0)) Chi << glob_av[ix]/(double)iblk << " " << err_x << " ";
  AveX.open("ave.chi."+to_string(metro),ios::app);
  if ((iblk == nblk) & (h == 0))  AveX << glob_av[ix]/(double)iblk << " " << err_x << " ";
  AveX.close();
  Chi.close();

  cout << "----------------------------" << endl << endl;
}

////////////////////////CONF. FINAL//////////////////////////////////

void ConfFinal(void){
  ofstream WriteConf;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)    WriteConf << s[i] << endl;
  WriteConf.close();
  rnd.SaveSeed();
}

////////////////////////PBC//////////////////////////////////
int Pbc(int i)  //Algorithm for periodic boundary conditions
{
  if(i >= nspin) i = i - nspin;
  else if(i < 0) i = i + nspin;
  return i;
}
////////////////////////ERROR//////////////////////////////////
double Error(double sum, double sum2, int iblk)
{
  if(iblk==1) return 0.0;
  else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
