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
double PSI (double,double,double);
double integral(double, double , double);


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

double x = 0; //starting position
double xnew  = 0; //new position
double p = 0, pp =0, q = 0, a=0;
int accepted = 0; //in order to have \sim 50% of acceptance ratio
int M = pow(10,6), L = 10000, N = int (M/L);
//int Neq = 1000;
double eps = 2.7; //step for markov chain metropolis
double h = 0;
double ave[N], ave2[N], sum[N], sum2[N], err[N]; //blocking average
ofstream print;
ofstream write;
ofstream stream;
print.open("print.out");
write.open("plot.out");
stream.open("position.out");

//ex :  --- fixed MU , SIGMA  (?)
//from the variational computation, we found out that the best parametres are:
//\mu = 0.815
//\sigma = 0.61

/////////////////////////////////////
//ex :  --- varying MU and SIGMA
double mu = 0.78, sigma = 0.58;
// if((mu == 0.78)&(sigma == 0.58))cout << "OK!" << endl;

for (int m = 1; m <= 20; m++){
	mu += 0.005;
	sigma = 0.58;
	for (int s = 1; s <= 20; s++){
		sigma += 0.005;
		accepted = 0;
		for (int i = 0; i < N; i++){
			ave[i] = 0;
			ave2[i] = 0;
			double accu = 0;
			for (int j=0; j<L; j++){
				xnew = rnd.Rannyu(x-eps,x+eps); //Uniform step, markov chain	
				p = PSI(mu, sigma, x);
				pp = PSI(mu, sigma, xnew);
				q = fmin( 1 , double(pp/p));
				a = rnd.Rannyu();
			
				if (a<q) { //accept
					x=xnew;
					accepted++;
					if (j%100 == 0) if (stream.is_open()) stream << x << " " ;
				} //else reject
				h = integral(mu,sigma,x);
				accu += h;
			}
			ave[i] = accu / L;
			ave2[i] = pow( ave[i] , 2);
	}
	//blocking uncertainties
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
	//printing
	cout << "mu = " << mu << " ; sigma = " << sigma << endl;
	if (print.is_open()) print << mu << " " << sigma << " " << sum[N-1] << " " << err[N-1] << " ";
	if (write.is_open()) if ((mu == 0.8 )&(sigma == 0.61+0.005))  for (int i=0; i<N; i++) write << sum[i] << " " << err[i]<< " " ;
	// if (write.is_open()) if ((mu == 0.815) & (sigma == 0.61))  write << sum[N-1] << " " << err[N-1]<< endl ; //for (int i=0; i<N; i++)
	cout << "Acceptance ratio = " << double(accepted)/M << endl;
	cout << "<H> = " << sum[N-1] << "; err = " << err[N-1] << endl;	
	}	
}
print.close();
write.close();
return 0;
}
///////////////////////////////////////////////////////

double error (double a, double b, int n){
	if (n==0) {return 0;}
	else {return sqrt((b-pow(a,2))/n);}
	}

double PSI (double  mu, double sigma, double x){
	double exp1 = exp(-pow(x-mu,2)/(2*pow(sigma,2)));
	double exp2 = exp(-pow(x+mu,2)/(2*pow(sigma,2)));      
	return pow(exp1+exp2,2);
}

double integral(double mu, double sigma, double x){
	double e1 = exp(-pow(x-mu,2)/(2*pow(sigma,2)));
	double e2 = exp(-pow(x+mu,2)/(2*pow(sigma,2))); 
    double n1 = pow(x-mu,2)*e1;
    double n2 = pow(x+mu,2)*e2;
    double d = e1+e2;
	return  pow(x,4) - 5./2.*pow(x,2)+1./(2*pow(sigma,2))-(n1+n2)/(2*pow(sigma,4)*d);
    // return -(n1+n2)/(d*(2*pow(sigma,2))) + pow(x,4) - 5./2.*pow(x,2);
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