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
   rnd.SaveSeed();
   
/////////////////////////////////////////////////// 
	//blocking
int M=pow(10,5);         // Total number of throws
int N=100;            // Number of blocks
int L=M/N; 	      // Number of throws in each block // L = 1000
double ave[N], ave2[N], sum[N], sum2[N], err[N];
double x;
	//data
double T = 1 ; 			//delivery time
double t = 0.01;		//time
double K = 100 ; 		//strike price
double r = 0.1 ; 		//risk-free interest rate
double mu = r;			//drift
double sigma = 0.25 ; 	//volatility
double S[100] ;			//price
double C = 0;			//Cost
double Ctrue[2] = {14.975790778311286 , 5.4595325819072364 }; //from Black-Scholes formulas:
S[0] = 100 ;			//asset price

int end = T/t ;

	//call & put directly
	for (int sign=0; sign < 2 ; sign++){ // sign = 0 --> call; sign = 1 --> put 
		for (int i=0; i<N; i++){
			ave[i] = 0;
			ave2[i] = 0;
			double accu = 0;
			for (int j=0; j<L; j++){
				x = rnd.Gauss(0,T);
				int index = T * 100 - 1;
				S[index] = S[0] * exp((mu - sigma*sigma / 2)*T + sigma * x );
				C = exp(-r*T) * max ( double(0) , pow(-1,sign) * S[index] + pow (-1,sign+1) * K ) ;
				accu += C;
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
		write.open("plot" + to_string(sign) + ".out");
		if (write.is_open())for (int i=0; i< N; i++)write << sum[i] - Ctrue[sign] << " " << err[i]<< " " ;
		else cerr << "PROBLEM: Unable to open random.out" << endl;
		write.close();
	}
	
	// call & put discretized

	for (int sign=2; sign < 4 ; sign++){ // sign = 0 --> call; sign = 1 --> put 
		for (int i=0; i<N; i++){
			ave[i] = 0;
			ave2[i] = 0;
			double accu = 0;
			for (int j=0; j<L; j++){
				for (int k=0; k < end-1 ; k++){
					x = rnd.Gauss(0,T);
					S[k+1] = S[k] * exp((mu - sigma*sigma / 2) * t + sigma * x *sqrt(t));
				}
				int index = T * T/t - 1;
				C = exp(-r*T) * max ( double(0) , pow(-1,sign) * S[index] + pow (-1,sign+1) * K ) ;
				accu += C;
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
		write.open("plot" + to_string(sign) + ".out");
		if (write.is_open())for (int i=0; i< N; i++)write << sum[i] - Ctrue[sign-2] << " " << err[i]<< " " ;
		else cerr << "PROBLEM: Unable to open random.out" << endl;
		write.close();
	}
	
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
