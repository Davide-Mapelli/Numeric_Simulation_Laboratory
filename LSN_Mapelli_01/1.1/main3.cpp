#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

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

   //for(int i=0; i<20; i++){ cout << rnd.Rannyu() << endl; }
   rnd.SaveSeed();

///////////part 3

int dataset[100];
double print[100];
double appo = 0;
double accu = 0;



	for (int j=0; j < 100; j++){
		accu = 0;
		for (int i = 0; i < 100; i++)dataset[i] = 0;
		
		for (int n = 0; n < pow(10,4); n++){
			appo = rnd.Rannyu();
			for (int i=0; i < 100; i++) {
				if ( (appo < (i+1)*1.f/100) & (appo > (i)*1.f/100) ){
					//cout << "appo = " << appo << endl ; 
					//cout << i*1.f/100 << endl;
					//cout << dataset[i] << endl;
					dataset[i]++;
					//cout << dataset[i] << endl;
				}
			}
		}
			
		for (int i = 0; i < 100; i++){
			accu += pow( dataset[i]-100 , 2) / 100; 
		}
		print[j] = accu ;
		
	}
	
   	ofstream write;
		write.open("plot3.out");
   		if (write.is_open()) for (int i=0; i< 100; i++) write << print[i] << " ";
   		else cerr << "PROBLEM: Unable to open random.out" << endl;
  	write.close();
	
return 0;
}