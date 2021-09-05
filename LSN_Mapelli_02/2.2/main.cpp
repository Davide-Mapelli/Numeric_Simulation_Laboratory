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
#include "posizione.h"

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
   
////////////////////////////////////////////////// 
int M=pow(10,4);      // Total number of throws
int N=100;           // Number of blocks
int L=M/N; 	     	 // Number of throws in each block // L = 1000

double ave[N], ave2[N], sum[N], sum2[N], err[N];
double print[100], printerr[100];	 //numers of steps on x-axis

////////part 1 - DISCRETE RW
int verso=0;
int dir=0;
	
	Posizione pos(0,0,0);
	
	for (int n=0; n<100; n++){
		print[n] = 0;
		printerr[n] = 0;
		//double accu = 0;
		for (int i=0; i<N; i++){
			double accu = 0;
			ave[i] = 0;
			ave2[i] = 0;
			for (int j=0; j<L; j++){
				pos.setX(0);
				pos.setY(0);
				pos.setZ(0);
				for (int m=0; m<n+1; m++){
					verso = rnd.Discrete(0,1)*2-1;
					dir = rnd.Discrete(0,2);
					if (dir == 0)pos.setX( pos.getX() + verso );
					if (dir == 1)pos.setY( pos.getY() + verso );
					if (dir == 2)pos.setZ( pos.getZ() + verso );
					//cout << pos.getR() << endl;
				}
				accu +=  pos.getR()*pos.getR();
			}
			ave[i] = accu / L;
			ave2[i] = pow(ave[i],2);
			//cout << i << " " << ave[i] << endl;
		}

		for (int i=0; i<N; i++){
			sum[i] = 0;
			sum2[i] = 0;
			//err[i] = 0;
			for(int j=0; j<i+1; j++){
				sum[i] += ave[j];
				sum2[i] += ave2[j];
			}
			sum[i] /= (i+1);
			sum2[i] /= (i+1);
			//cout << i << " " << sum[i] << endl;
			err[i] = error (sum[i], sum2[i], i);
			
		}
		print[n] = sum[N-1];
		//cout << print [n] << endl;
		printerr[n] = err[N-1];
	}

	ofstream write;
	write.open("plot1.out");
   	if (write.is_open())for (int i=0; i< 100; i++)write << sqrt(print[i]) << " " << sqrt(printerr[i]) << " " ;
   	else cerr << "PROBLEM: Unable to open random.out" << endl;
  	write.close();
	
////////part 2 - CONTINUE RW
double theta=0; // [0, pi]
double phi = 0; // [0, 2pi]

Posizione pos1(0,0,0);
	
	for (int n=0; n<100; n++){
		print[n] = 0;
		printerr[n] = 0;
		//double accu = 0;
		for (int i=0; i<N; i++){
			double accu = 0;
			ave[i] = 0;
			ave2[i] = 0;
			for (int j=0; j<L; j++){
				pos1.setX(0);
				pos1.setY(0);
				pos1.setZ(0);
				for (int m=0; m<n+1; m++){
					theta = acos(1 - 2*rnd.Rannyu());
					phi = rnd.Rannyu(0, 2*M_PI);
					pos1.setX( pos1.getX() + sin(theta) * cos(phi));
					pos1.setY( pos1.getY() + sin(theta) * sin(phi));
					pos1.setZ( pos1.getZ() + cos(theta));
				}
				accu +=  pos1.getR()*pos1.getR();
			}
			ave[i] = accu / L;
			ave2[i] = pow(ave[i],2);
			//cout << i << " " << ave[i] << endl;
		}

		for (int i=0; i<N; i++){
			sum[i] = 0;
			sum2[i] = 0;
			//err[i] = 0;
			for(int j=0; j<i+1; j++){
				sum[i] += ave[j];
				sum2[i] += ave2[j];
			}
			sum[i] /= (i+1);
			sum2[i] /= (i+1);
			//cout << i << " " << sum[i] << endl;
			err[i] = error (sum[i], sum2[i], i);
			
		}
		print[n] = sum[N-1];
		//cout << print [n] << endl;
		printerr[n] = err[N-1];
	}

	ofstream write1;
	write1.open("plot2.out");
   	if (write1.is_open())for (int i=0; i< 100; i++)write1 << sqrt(print[i]) << " " << sqrt(printerr[i]) << " " ;
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
