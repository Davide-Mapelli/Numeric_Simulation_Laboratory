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
   
  /////////////////////////////////////////////////// BUFFON
   
	double d = 1;
	double l = 0.5;

	int M=1000000;         // Total number of throws
	int N=100;            // Number of blocks
	int L=M/N; 	          // Number of throws in each block // L = 1000

	int hit = 0;
	double a=0;
	double b=0;
	double pi = 0;
	
	double ave[N];
	double ave2[N];
	double sum[N];
	double sum2[N];
	double err[N];

	double x,y;
	//double accu = 0;

	for (int i=0; i<N; i++){
		hit=0;
		ave[i] = 0;
		ave2[i] = 0;
		for (int j=0; j<L; j++){
			a = rnd.Rannyu(0,10);
			do {
				x = rnd.Rannyu(a-l,a+l);
				y = rnd.Rannyu(-l,l);
				}
			while((sqrt(pow(x-a,2)+pow(y,2)) > l));
			b = a + l * (x-a)/sqrt(pow(x-a,2)+pow(y,2));
			if (floor(a)!= floor(b))hit +=1;
		}
	pi = (2*l*L) /(hit *d);
	ave[i] = pi;
	ave2[i] = pow(pi , 2);
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

		for (int i=0; i< N; i++){cout << sum[i] -  M_PI << " " << err[i]<< endl ;}
 		
  		ofstream write;
		write.open("plot.out");
   		if (write.is_open()){
   			for (int i=0; i< N; i++){write << sum[i] - M_PI << " " << err[i]<< " " ;}
   					}
   		else cerr << "PROBLEM: Unable to open random.out" << endl;
  		write.close();
   				
   
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
