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

using namespace std;
 
double error (double a, double b, int n){
	if (n==0) {return 0;}
	else {return sqrt((b-pow(a,2))/n);}
	}

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

   //for(int i=0; i<20; i++){ cout << rnd.Rannyu() << endl; }
   rnd.SaveSeed();
   
 /////////////////////////////////////////////////// 
  
int M=100000;         // Total number of throws
int N=100;            // Number of blocks
int L=M/N; 	      // Number of throws in each block // L = 1000
////////part 1
double ave[N], ave2[N], sum[N], sum2[N], err[N];
double x;

	for (int i=0; i<N; i++){
		ave[i] = 0;
		ave2[i] = 0;
		double accu = 0;
		for (int j=0; j<L; j++){
			x = rnd.Rannyu();
			accu += M_PI/2 * cos(M_PI*x/2);
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
	
	ofstream write;
	write.open("plot1.out");
   	if (write.is_open())for (int i=0; i< N; i++)write << sum[i] - 1 << " " << err[i]<< " " ;
   	else cerr << "PROBLEM: Unable to open random.out" << endl;
  	write.close();
	
////////part 2//Importance sampling


double y;
	for (int i=0; i<N; i++){
		ave[i] = 0;
		ave2[i] = 0;
		double accu = 0;
		for (int j=0; j<L; j++){
			x = rnd.Rannyu();
			y = 1 - sqrt(1-x);
			accu += M_PI/2 * cos(M_PI*y/2) / (2*(1-y));
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

	ofstream write1;
	write1.open("plot2.out");
   	if (write1.is_open())for (int i=0; i< N; i++)write1 << sum[i] - 1 << " " << err[i]<< " " ;
   	else cerr << "PROBLEM: Unable to open random.out" << endl;
  	write1.close();

return 0;
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
