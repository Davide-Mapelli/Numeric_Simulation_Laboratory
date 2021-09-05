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

//Cities creation:
unsigned int cs = 32; //number of cities
city Sim(cs, rnd);
unsigned int N = 1; //member of the population
bool type = 0; //0 --> CIRCLE; 1 --> SQUARE
for (int c =0; c<2; c++ ){
	 
	// unsigned int c = 32; //number of cities
	if(type == 0)Sim.circle1D();
	if(type == 1)Sim.square2D();
	//printing cities
	for (unsigned int i=0; i<Sim.GetN(); i++) cout << "città n° " << i+1 << ": ( " << Sim.GetPosX(i) << " ; " << Sim.GetPosY(i) <<  " )" << endl;
	Sim.print();
	//generating pop
	GAs pop(Sim.GetN(), N, Sim, rnd); //every row must have m_N number of column, first constraint. I insert 100 populator in my pop ( in reality they should be NxN)
	pop.Check();
	cout << endl;	
	//printing pop
	for (unsigned int i = 0; i<N; i++){
		for (unsigned int j = 0; j<Sim.GetN(); j++){
			cout << pop.GetPos(pop.GetMember(i),j) << ", " ;
		}
		cout << endl;
	}
	cout << pop.Check() << endl;

	for (unsigned int i=0; i<N; i++) {
		for	(unsigned int j=0; j<Sim.GetN(); j++){
			cout << pop.GetPos(pop.GetMember(i),j) << " , " ;		
		}
		cout << endl;
		cout << "L_1 = " << pop.NormL(pop.GetMember(i), 0);
		cout << " ; L_2 = " << pop.NormL(pop.GetMember(i), 1) << endl;
	}

	//sorting paths	
	pop.Sorter();
	
	cout << "... sorting chromosomes ... " << endl;
	
	for (unsigned int i=0; i<pop.GetR(); i++) {
		for	(unsigned int j=0; j<pop.GetC(); j++){
			cout << pop.GetPos(pop.GetMember(i),j) << " , " ;		
		}
		cout << endl;
		cout << "L_1 = " << pop.NormL(pop.GetMember(i), 0);
		cout << " ; L_2 = " << pop.NormL(pop.GetMember(i), 1) << endl;
	}
	
	cout << pop.Check() << endl;

///////////////////////////////////////
//SIMULATED ANNEALING
///////////////////////////////////////	
double beta = 0.01;
double temp = 1./beta;
double deltabeta = 0.001;
double betamax = 50; 

double deltaL = 0, Lold =0, Lnew = 0;
double p=0,q=0;
int n;

pop.SetPath(pop.GetMember(N-1));

Lold = pop.NormL(pop.GetPath(),0);
ofstream print;
print.open("Length" + to_string(type) + ".out");	
while (beta < betamax){
	if (print.is_open())print << Lold << " ";

	for (int k = 0; k < 256; k++){
	//create a new path
	n = rnd.Discrete(0,3);
	if (n==0)pop.PairPermutation();
	if (n==1)pop.Shift();
	if (n==2)pop.Permutation();
	if (n==3)pop.Inversion();

	Lnew = pop.NormL(pop.GetNewPath(),0);
	deltaL =  Lnew - Lold;

	if (deltaL<0){
		Lold = Lnew;
		pop.SetPath(pop.GetNewPath());
	}
	else{
		p = exp(-beta*deltaL);
		q = rnd.Rannyu();
		if (q<p){
			Lold = Lnew;
			pop.SetPath(pop.GetNewPath());
		}
	}
	}
	beta += deltabeta;
	temp = 1./beta;

	cout << "beta = " << beta << "; temp = "<< temp << endl;
	// for (unsigned int i=0; i < cs; i++)cout << pop.GetPath()[i] << " " ;

}

ofstream write;
	if (type==0)write.open("citycircle.out");	
	if (type==1)write.open("citysquare.out");
	for (unsigned int k=0; k<cs; k++){
		//for (unsigned int i=0; i<Sim.GetN(); i++) 
		unsigned int i = pop.GetPath()[k];
		cout << "i = " << i;
		cout << " , città n° " << k+1 << ": ( " << Sim.GetPosX(i-1) << " ; "<< Sim.GetPosY(i-1) <<  " )" << endl;
		if (write.is_open())write << Sim.GetPosX(i-1) << " " << Sim.GetPosY(i-1) << " " ;
		else cerr << "PROBLEM: Unable to open random.out" << endl;
	}
	write.close();
	print.close();

	type += 1;
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
