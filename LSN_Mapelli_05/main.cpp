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
#include <string>
#include <cmath>
#include "random.h"

double error (double, double, int);
double PSI (double [3], int );

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   rnd.SaveSeed();
   
/////////////////////////////////////////////////// 
//BLOCKING

double position[3]; //starting position
double xyz[3]; //new position
double p = 0, pp =0, q = 0, x=0;
int accepted = 0;
int M = pow(10,6), L = 10000, N = int (M/L), Neq = 100;
double epsilon[2] = {1.2, 2.9};
double sigma[2] = {0.75, 1.85};
double multivariate[6] = {1.25, 0.75, 0.25, 2.85,1.85,0.85};


double Rstart[2] = {0,100};
double r, R;
double ave[N], ave2[N], sum[N], sum2[N], err[N];
string step[3]= {"Uniform", "Gauss", "multivariate"};
ofstream print;
ofstream stream;

for (int rr = 0; rr < 2; rr++){
	R = Rstart[rr]; //r start
	cout << "R_start = " << R << " a_0" <<  endl;

	for (int s = 0; s < 3; s++){

		cout << "Step : " << step[s] << endl;

		for (int n = 1; n <= 2 ; n++){ // n = 1; n = 2
		
			for (int k=0; k<3; k++)position[k] = sqrt(double(R*R)/3); //initialization
			accepted = 0;
			if (R == 0) print.open("print." + step[s] + to_string(n) + ".out"); //grafico 3D
			
			/////////////////
			//EQUILIBRATION//
			/////////////////

			cout << "Equilibration phase: n = " << n << endl;
			//if (Neq > 0){cout << "Equilibration phase: n = " << n << endl;}
			for (int i = 0; i < Neq; i++){
				if (s == 0) for (int k=0; k<3; k++)xyz[k] = rnd.Rannyu(position[k]-epsilon[n-1],position[k]+epsilon[n-1]); //Uniform step	
				if (s == 1) for (int k=0; k<3; k++)xyz[k] = rnd.Gauss(position[k],sigma[n-1]); //Gauss step
				if (R == 0) if (s == 2) for (int k=0; k<3; k++)xyz[k] = rnd.Gauss(position[k],multivariate[(n-1)*3+k]);//Multivariate normal transition probability (Gaussian for each coordinate) step [only origin as starting point]
				p = PSI(position,n);
				pp = PSI(xyz,n);
				q = fmin( 1, double(pp/p));
				x = rnd.Rannyu();
				if (x<q) { //accept
					for (int k=0; k<3; k++)position[k]=xyz[k];
					//accepted++;
					} 
				//reject
				}

			///////////////
			//MEASUREMENT//
			///////////////

			//cout << "Measurement phase: n = " << n << endl;
			for (int i = 0; i < N; i++){
				ave[i] = 0;
				ave2[i] = 0;
				double accu = 0;
				for (int j=0; j<L; j++){
					if (s == 0) for (int k=0; k<3; k++)xyz[k] = rnd.Rannyu(position[k]-epsilon[n-1],position[k]+epsilon[n-1]); //Uniform step	
					if (s == 1) for (int k=0; k<3; k++)xyz[k] = rnd.Gauss(position[k],sigma[n-1]); //Gauss step	
					if (R == 0) if (s == 2) for (int k=0; k<3; k++)xyz[k] = rnd.Gauss(position[k],multivariate[(n-1)*3+k]);
					p = PSI(position,n);
					pp = PSI(xyz,n);
					q = fmin( 1 , double(pp/p));
					x = rnd.Rannyu();
					
					if (x<q) { //accept
						for (int k=0; k<3; k++)position[k]=xyz[k];
						accepted++;
					} //else reject
					if (j%100 == 0)if (print.is_open())for (int k=0; k< 3; k++)print << position[k] << " " ;
					r = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
					accu += r;
				}
				ave[i] = accu / L;
				ave2[i] = pow( ave[i] , 2);
			}
			
			for (int i=0; i<N; i++){
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
			
			print.close();

			ofstream write;
			if (R == 0) write.open("plot." + step[s] + to_string(n) + ".out");
			if (R != 0) write.open("plot.far" + step[s] + to_string(n) + ".out");
			if (write.is_open())for (int i=0; i< N; i++)write << sum[i] << " " << err[i]<< " " ;
			write.close();
			cout << "Acceptance ratio = " << double(accepted)/M << endl;
		}
	}	
}
	return 0;
}
///////////////////////////////////////////////////////

double error (double a, double b, int n){
	if (n==0) {return 0;}
	else {return sqrt((b-pow(a,2))/n);}
	}

double PSI (double x[3], int n){ //modulo quadro della funzione d'onda
	double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	if (n==1)return exp(-2*r)/M_PI; //100
	else {return (1/(32*M_PI)) * exp(-r) * x[2]*x[2];  } //210 //cos^2	theta * r^2 = z^2
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