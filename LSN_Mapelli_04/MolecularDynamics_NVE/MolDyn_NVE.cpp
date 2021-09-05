/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include "MolDyn_NVE.h"

using namespace std;

/*
  bool type;
  cout << "Type 0 if you want to start from a new config (config.0)" << endl; 
  cout << "Type 1 if you want to start from an old config (old.final)" << endl;
  cin >> type ;
*/

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////MAIN////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

int main(){ 
  Input();  //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
    Move();           //Move particles with Verlet algorithm
    if(istep%iprint == 0) cout << "Number of time-steps: " << istep << " / " << nstep << endl;
    if(istep%10 == 0){
      Measure();           //Properties measurement
      //  ConfXYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
      nconf += 1;
    }
   if (istep==nstep-1)ConfFinal(1); //old.0
  }
  ConfFinal();           //Write final configuration to restart : config.final --> old.final
  Average();
  //Gofr();
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////FUNCTIONS////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

//////////////INPUT////////////////////////////////

void Input(){ //Prepare all stuff for the simulation

  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> type;
  cout << "Type: " << type << endl;
  ReadInput >> temp;
  ReadInput >> tempstar;
  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;
  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie =  2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  //measurement of g(r) : g of r
  igofr = n_props;
  // nbins = 100;
  n_props = igofr + nbins;
  for (int i =0; i < nbins; i++) bin[i] = 0;
  bin_size = (box*sqrt(3)/2.0)/(double)nbins;

/////NEW CONFIG/////
  if (type == 0){ 

   //Read initial configuration//
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
   }
    ReadConf.close();
  
    //Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }

   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }

   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }

   //return;
  }

/////OLD CONFIG//////

 else{
  //Read initial configuration
   cout << "Read initial configuration from file old.final (𝑟⃗ (𝑡)) " << endl << endl;
   ReadConf.open("old.final");
   for (int i=0; i<npart; ++i){
     ReadConf >> x[i] >> y[i] >> z[i];
     x[i] = x[i] * box;
     y[i] = y[i] * box;
     z[i] = z[i] * box;
   }
   ReadConf.close(); 

   cout << "Read previous configuration from file old.0 (𝑟⃗ (𝑡−𝑑𝑡)) " << endl << endl;
   ReadConf.open("old.0");
   for (int i=0; i<npart; ++i){
     ReadConf >> xold[i] >> yold[i] >> zold[i];
     xold[i] = xold[i] * box;
     yold[i] = yold[i] * box;
     zold[i] = zold[i] * box;
   }
   ReadConf.close();

//Prepare initial velocities
   cout << "Prepare velocities using 𝑟⃗(𝑡−𝑑𝑡) and 𝑟⃗(𝑡) ... " << endl << endl;
   Move();

   for (int i=0; i<npart; ++i){
     vx[i] = (Pbc(x[i] - xold[i]))/(0.5 * delta);
     vy[i] = (Pbc(y[i] - yold[i]))/(0.5 * delta);
     vz[i] = (Pbc(z[i] - zold[i]))/(0.5 * delta);
   }    
  
  //Rescale velocities
  if (tempstar!=0){
    cout << "Rescaling velocities using 𝑟⃗(𝑡+𝑑𝑡) and 𝑟⃗(𝑡) and temperature " << endl << endl;
    double K, Teff;
    for (int i=0; i<npart; ++i) K += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    Teff = (2.0 / 3.0) * K/(double)npart; //Temperature

    for (int i=0; i<npart; ++i){
     vx[i]*= tempstar/Teff;
     vy[i]*= tempstar/Teff;
     vz[i]*= tempstar/Teff;
     }    
    for (int i=0; i<npart; ++i){
     xold[i] = Pbc ( x[i] - delta * vx[i] ) ;
     yold[i] = Pbc ( y[i] - delta * vy[i] ) ;
     zold[i] = Pbc ( z[i] - delta * vz[i] ) ;
     }
  }
 }

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Measure(){ //Properties measurement
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
  string state;

  if (rho == 1.1)state = "solid";
  if (rho == 0.8)state = "liquid";
  if (rho == 0.05)state = "gas";

  
  Epot.open("output_"+state+"_epot.dat",ios::app);
  Ekin.open("output_"+state+"_ekin.dat",ios::app);
  Temp.open("output_"+state+"_temp.dat",ios::app);
  Etot.open("output_"+state+"_etot.dat",ios::app);
  
  
  v = 0.0; //reset observables
  t = 0.0;

  //reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) stima_gofr[k]=0.0;


  //  cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

    
    // g(r)
    if (dr <= box/2){
      for(int k=0; k<nbins; k++){
        if ( k == int(dr/bin_size) )stima_gofr[igofr+k] += 2;
      }
    } 

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

 //Potential energy
       v += vij;
     }
    }          
  }
 //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  
  
  //writing instant values 
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  

  for(int i=0; i<nbins; ++i)blk_av[i] += bin[i];
  blk_norm += 1.0;
  norm += 1.0;
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  
cout << "Print final configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfFinal(bool type){ //Write final configuration
  
  ofstream WriteConf;

  cout << "Print final configuration to file old.0" << endl;
  WriteConf.open("old.0");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Block Averaging

void Average (void){
  //int M=nstep;           // Total number of throws
  int L=20;           //Length of the block
  //int N=M/L;            // Number of blocks
  double x;
  int M = 0;
  int conta = 0;
  string title[4] = {"ekin", "epot", "etot", "temp"};
  //double ave[4][N], ave2[4][N], sum[4][N], sum2[4][N], err[4][N];
  string state;
  if (rho == 1.1)state = "solid";
  if (rho == 0.8)state = "liquid";
  if (rho == 0.05)state = "gas";

  ifstream Read;
  Read.open("output_"+state+"_ekin.dat"); //Read input
  while(!Read.eof()){
    Read >> x;
    M++;
    //cout << N;
  }
  Read.close();
  int N = M/L;
  double ave[N], ave2[N], sum[N], sum2[N], err[N];
  for (int property = 0; property < 4 ; property++){ // EKin, Epot, Etot, Temp
    ifstream ReadData;
    ReadData.open("output_"+ state +"_"+ title[property]+ ".dat"); //Read input
    
		for (int i = 0; i<50; i++){
			ave[i] = 0;
			ave2[i] = 0;
			double accu = 0;
			for (int j=0; j<L; j++){ //L=20
        while(conta < (M-1000)){
          ReadData >> x;
          conta++;
        }
				ReadData >> x;
				accu += x;
			}
			ave[i] = accu / L;
			ave2[i] = pow( ave[i] , 2);
      //i++;
		}
		for (int i=0; i<50; i++){
			sum[i] = 0;
			sum2[i] = 0;
			err[i] = 0;
			for(int j=0; j<i+1; j++){
				sum[i] += ave[j];
				sum2[i] += ave2[j];
				}
			sum[i] /= (i+1);
			sum2[i] /= (i+1);
			err[i] = error (sum[i], sum2[i], i);
		}

		ofstream write;
		write.open("ave_" + state +"_"+ title[property] + ".out");
		if (write.is_open())for (int i=0; i< 50; i++)write << sum[i] << endl << err[i]<< endl;
		else cerr << "PROBLEM: Unable to open random.out" << endl;
		write.close();
    ReadData.close();   

    ofstream Gave;
    Gave.open("output." + state + ".gave.0",ios::app);
    for(int i=igofr; i<n_props;i++){
      deltaV = 4./3.* M_PI * ( pow((i-igofr+1)*bin_size,3) - pow((i-igofr)*bin_size,3) );
      stimagofr = blk_av[i]/blk_norm / (rho * m_part * deltaV);
      Gave << "\t" << stimagofr;
      glob_av[i] += stimagofr;
      glob_av2[i] += stimagofr*stimagofr;
    }
    Gave << endl;
    Gave.close();
	}
}
//////////////////////////////
/*
void Gofr (void){

  ofstream Gofr;
  Gofr.open("output.gofr.0");
  double coeff = ((double)npart*rho*vol);
  //for (int i=0; i < nbins; i++)Gofr << bin[i]/coeff  << " ";

  ofstream Gave;
  Gave.open("output.gave.0");
  for (int i=0; i < nbins; i++){
    stima_gofr[i] = blk_av[i]/coeff;
    glob_av[i] += stima_gofr[i];
    glob_av2[i] += stima_gofr[i]*stima_gofr[i];
    err_gofr[i]=error(glob_av[i],glob_av2[i],iblk);
    Gofr << glob_av[i]/(double)iblk  << " ";
    Gave << glob_av[i]/(double)iblk << " " << err_gofr[i] << " "; //if (iblk == nblk) 
  }
  Gofr.close();
  Gave.close();
}
*/

/*
//g(r)
    
    for(int i=igofr; i<n_props;i++){
      deltaV = 4./3.* pi * ( pow((i-igofr+1)*bin_size,3) - pow((i-igofr)*bin_size,3) );
      stima_gofr = blk_av[i]/blk_norm / (rho * m_part * deltaV);
      Gofr << "\t" << stima_gofr;
      glob_av[i] += stima_gofr;
      glob_av2[i] += stima_gofr*stima_gofr;
    }
    Gofr << endl;
    if (iblk==nblock){
		for (int i = igofr; i < n_props; i++){
      		err_gofr = Error(glob_av[i], glob_av2[i], nblock);
      		Gave << "\t" << (i-igofr+1./2.)*bin_size << "\t" << glob_av[i] / nblock << setw(wd) << err_gofr << endl;
		}
	}



*/
////////////////////////////////////////////////////////////

double error (double a, double b, int n){
	if (n==0) {return 0;}
	else {return sqrt((b-pow(a,2))/n);}
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
