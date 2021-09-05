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
double accu = 0;
int reps[4]= {1,2,10,100};
int N=0;

ofstream write;
write.open("plot1.out");
	for (int n=0; n<4; n++){
		for (int j=0; j< pow(10,4); j++){
			N = reps[n];
			accu = 0;
			for (int i= 0; i<N; i++)accu += rnd.Rannyu();
			if (write.is_open()) write <<  accu * 1.f / N << " " ;
			else cerr << "PROBLEM: Unable to open random.out" << endl;
		}
	}
write.close();

ofstream write1;
write1.open("plot2.out");
	for (int n=0; n<4; n++){
		for (int j=0; j< pow(10,4); j++){
			N = reps[n];
			accu = 0;
			for (int i= 0; i<N; i++)accu += rnd.Exponential(1);
			if (write1.is_open()) write1 <<  accu * 1.f / N << " " ;
			else cerr << "PROBLEM: Unable to open random.out" << endl;
		}
	}
write1.close();

ofstream write2;
write2.open("plot3.out");
	for (int n=0; n<4; n++){
		for (int j=0; j< pow(10,4); j++){
			N = reps[n];
			accu = 0;
			for (int i=0; i<N; i++) accu += rnd.Cauchy_Lorentz(0,1);
			if (write2.is_open()) write2 <<  accu * 1.f / N << " " ;
			else cerr << "PROBLEM: Unable to open random.out" << endl;
		}
	}
write2.close();

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
