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
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include "random.h"
#include "GAs.h"
#include "city.h"

using namespace std;
 
 bool myfunction (GAs M1, unsigned int i, unsigned int j) {
	vector <unsigned int> m1 = M1.GetMember(i);
	vector <unsigned int> m2 = M1.GetMember(j);
	return  M1.NormL(m1,0) < M1.NormL(m2,0);
 }

class mycomparator {
public : 
	bool operator() (GAs M1, vector <unsigned int> m1, vector <unsigned int> m2){
		return M1.Comparator (m1, m2) ;
		};
} funct ;

int main (int argc, char *argv[]){
	//INITIALIZATION PSEUDO-CASUAL NUMBER GENERATOR CLASS
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
/////////////////CITIES///////////////////

unsigned int cs = 32; //number of cities
city Sim(cs, rnd);
unsigned int type = 0; // 0--> Cirle; 1 --> Square
unsigned int N = 256; //member of the population

for(int r = 0 ; r < 2 ; r++){
	
	if(type == 0) Sim.circle1D();
	if(type == 1) Sim.square2D();

	ofstream print;
	print.open("Length"+to_string(type)+".out");

	//print cities
	//for (unsigned int i=0; i<cs; i++) cout << "città n° " << i+1 << ": ( " << Sim.GetPosX(i) << " ; " << Sim.GetPosY(i) <<  " )" << endl;
	Sim.print();
	GAs population(cs, N, Sim, rnd); //every row must have m_N number of column, first constraint. I insert 100 populator in my population ( in reality they should be NxN)
	population.Check();
	cout << endl;
	//print chromosomes
	/*
	for (unsigned int i = 0; i<N; i++){
		for (unsigned int j = 0; j<Sim.GetN(); j++)cout << population.GetPos(population.GetMember(i),j) << ", " ;
		cout << endl;
	}
	cout << population.Check() << endl;
	for (unsigned int i=0; i<N; i++) {
		for	(unsigned int j=0; j<cs; j++)cout << population.GetPos(population.GetMember(i),j) << " , " ;
		cout << endl;
		cout << "L_1 = " << population.NormL(population.GetMember(i), 0);
		cout << " ; L_2 = " << population.NormL(population.GetMember(i), 1) << endl;
	}
	*/

	//print sorted chromosomes
	cout << "... sorting chromosomes ... " << endl;
	/*
	for (unsigned int i=0; i<population.GetR(); i++) {
		for	(unsigned int j=0; j<population.GetC(); j++){
			cout << population.GetPos(population.GetMember(i),j) << " , " ;		
		}
		cout << endl;
		cout << "L_1 = " << population.NormL(population.GetMember(i), 0);
		cout << " ; L_2 = " << population.NormL(population.GetMember(i), 1) << endl;
	}
	*/
	//cout << population.Check() << endl;
	//cout << "selected : " << population.Selection() << endl;
	//for (int i=0; i< 100; i++)population.Permutation();
	/*
	for (unsigned int i=0; i<N; i++) {
		for	(unsigned int j=0; j<cs; j++){
			cout << population.GetPos(population.GetMember(i),j) << " , " ;		
		}
		cout << endl;
	}
	*/
	population.Sorter();
	/*
	for (unsigned int i=0; i<N; i++) {
		cout << "L_1 = " << population.NormL(population.GetMember(i), 0);
		cout << " ; L_2 = " << population.NormL(population.GetMember(i), 1) << endl;
	}
	*/
	cout << "-----------" << endl;
	

//////////////////////////////////
///////GENETIC ALGORITHMS/////////
//////////////////////////////////
	double x = 0;
	unsigned int n = 0;
	vector <vector<unsigned int>> appo;
	unsigned int gen0 = 512;
	unsigned int gen1 = 1024;

	unsigned int gen = 0;
	if (type == 0)gen = gen0;
	if (type == 1)gen =gen1;

	for (unsigned int j=0; j<gen; j++) { //n° of generations
		unsigned int i = 0;
		while (i < N){
			if (i != N-1){
			x = rnd.Rannyu(0,1);
			if(x<0.4){
				n = rnd.Discrete(0,3);
				if (j==128)cout << " n  = " << n << endl;
				if (n==0)population.PairPermutation(i);
				if (n==1)population.Shift(i);
				if (n==2)population.Permutation(i);
				if (n==3)population.Inversion(i);
				i += 1;
			}
			else{
				population.Crossover(i,i+1);
				i += 2;
			}
			}
			else{
				x = rnd.Rannyu(0,1);
				if(x<0.4){
					n = rnd.Discrete(0,3);
					if (n==0)population.PairPermutation(i);
					if (n==1)population.Shift(i);
					if (n==2)population.Permutation(i);
					if (n==3)population.Inversion(i);
					i+=1;
				}
			}
		}		

		if (j%16 == 0 ) cout << "n = " << j << endl;
		population.SetPopulation(population.GetSon());
		population.Sorter();

		double accu = 0;
		for(unsigned int p = N/2; p < N; p++)accu += population.NormL(population.GetMember(p), 0);
		if (print.is_open()) print << 2.*accu/N << " " << population.NormL(population.GetMember(N-1), 0) << " " ;
	}
/*
	for (unsigned int i=0; i<N; i++) {
		cout << "L_1 = " << population.NormL(population.GetMember(i), 0);
		cout << " ; L_2 = " << population.NormL(population.GetMember(i), 1) << endl;
	}
	*/
	ofstream write;
	if(type == 0)write.open("citycircle.out");	
	if(type == 1)write.open("citysquare.out");
	for (unsigned int k=0; k<cs; k++){
		//for (unsigned int i=0; i<Sim.GetN(); i++) 
		unsigned int i = population.GetPos(population.GetMember(N-1),k);
		cout << "i = " << i;
		cout << " , città n° " << k+1 << ": ( " << Sim.GetPosX(i-1) << " ; "<< Sim.GetPosY(i-1) <<  " )" << endl;
		if (write.is_open())write << Sim.GetPosX(i-1) << " " << Sim.GetPosY(i-1) << " " ;
		else cerr << "PROBLEM: Unable to open random.out" << endl;
	}
	write.close();
	print.close();
	type +=1 ;
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
