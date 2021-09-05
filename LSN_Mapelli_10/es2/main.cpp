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
#include <vector>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <iomanip>
#include "mpi.h"
#include "random.h"
#include "GAs.h"
#include "city.h"

using namespace std;

/*
MPI_Init() -> MPI initialization
MPI_Comm_size() -> to know how many processes
MPI_Comm_rank() -> to know who am I
MPI_Send() -> to send a message
MPI_Recv() -> to receive a message
MPI_Finalize() -> MPI finalization
*/
template <class T>
T *resize(vector<vector<T>> mat)
{
	T *arr = new T[mat.size() * mat[0].size()];
	for (int i = 0; i < mat.size(); i++)
		for (int j = 0; j < mat[i].size(); j++)
			arr[i * mat[0].size() + j] = mat[i][j];
	return arr;
};
template <class T>
vector<vector<T>> resize(T *v, int size, int rowsize)
{
	vector<vector<T>> mat(int(size / rowsize), vector<T>(rowsize, 0));
	for (int i = 0; i < mat.size(); i++)
		for (int j = 0; j < mat[i].size(); j++)
			mat[i][j] = v[i * rowsize + j];
	return mat;
};

 
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

	/////// INITIALIZATION MPI
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat;
	MPI_Request req;

	//INITIALIZATION PSEUDO-CASUAL NUMBER GENERATOR CLASS
	Random rnd("seed" + to_string(rank) + ".in", "seed" + to_string(rank) + ".out");
	
unsigned int cs = 32; //number of cities
unsigned int N = 256; //member of the pop
unsigned int gen = 512; //number of generations
// double *map1D;
// int *world1D;
// map1D = new double[2*cs]; //x,y
// world1D = new int[cs*N*size];
//vector<vector<int>> world;//
city Sim(cs, rnd);
//leggi città da testo 
ifstream incities;
incities.open("city.out");

for (unsigned int i = 0; i < cs; i++){
	double appox=0;
	double appoy=0;
	incities >>  appox  >> appoy;
	Sim.SetPosX(i,appox);
	Sim.SetPosY(i,appoy);
}
incities.close();
	// MPI_Bcast(map1D, 2 * cs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// if (rank != 0) map = resize<double>(map1D, 2 * cs, 2);
	GAs pop(cs, N, Sim, rnd); //every row must have m_N number of column, first constraint. I insert 100 populator in my pop ( in reality they should be NxN)
	pop.Check();
	/*
	int *pop1D;
	pop1D = new int[N*cs];
	for (unsigned int q = 0; q<N*cs;q++)pop1D[q] = pop.GetMember(q/cs)[q%cs];
	
	MPI_Gather(pop1D, N * cs, MPI_INTEGER, world1D, N * cs, MPI_INTEGER, 0, MPI_COMM_WORLD);
	delete[] pop1D;
	*/
	
	// if (rank == 0)pop.Sorter();
		//world.assign(world.begin(), world.begin() + nroutes);
	
	//print cities
	for (unsigned int i=0; i<cs; i++) cout << "città n° " << i+1 << ": ( " << Sim.GetPosX(i) << " ; " << Sim.GetPosY(i) <<  " )" << endl;
	//Sim.print();
	
	cout << endl;
	//print chromosomes
	/*
	for (unsigned int i = 0; i<N; i++){
		for (unsigned int j = 0; j<Sim.GetN(); j++)cout << pop.GetPos(pop.GetMember(i),j) << ", " ;
		cout << endl;
	}
	cout << pop.Check() << endl;
	for (unsigned int i=0; i<N; i++) {
		for	(unsigned int j=0; j<cs; j++)cout << pop.GetPos(pop.GetMember(i),j) << " , " ;
		cout << endl;
		cout << "L_1 = " << pop.NormL(pop.GetMember(i), 0);
		cout << " ; L_2 = " << pop.NormL(pop.GetMember(i), 1) << endl;
	}
	*/

	//print sorted chromosomes
	cout << "... sorting chromosomes ... " << endl;
	/*
	for (unsigned int i=0; i<pop.GetR(); i++) {
		for	(unsigned int j=0; j<pop.GetC(); j++){
			cout << pop.GetPos(pop.GetMember(i),j) << " , " ;		
		}
		cout << endl;
		cout << "L_1 = " << pop.NormL(pop.GetMember(i), 0);
		cout << " ; L_2 = " << pop.NormL(pop.GetMember(i), 1) << endl;
	}
	*/
	//cout << pop.Check() << endl;
	//cout << "selected : " << pop.Selection() << endl;
	//for (int i=0; i< 100; i++)pop.Permutation();
	/*
	for (unsigned int i=0; i<N; i++) {
		for	(unsigned int j=0; j<cs; j++){
			cout << pop.GetPos(pop.GetMember(i),j) << " , " ;		
		}
		cout << endl;
	}
	*/
	pop.Sorter();
	/*
	for (unsigned int i=0; i<N; i++) {
		cout << "L_1 = " << pop.NormL(pop.GetMember(i), 0);
		cout << " ; L_2 = " << pop.NormL(pop.GetMember(i), 1) << endl;
	}
	*/
	// cout << "-----------" << endl;
	

//////////////////////////////////
///////GENETIC ALGORITHMS/////////
//////////////////////////////////
	double x = 0;
	unsigned int n = 0;
	vector <vector<unsigned int>> appo;
	vector <unsigned int> appo2;
	unsigned int *migr_send = new unsigned int[cs * 3];
	unsigned int *migr_rec = new unsigned int[cs * 3];
	int shift = 1;
	
	ofstream print;
	print.open("Length"+to_string(rank)+".out");

	for (unsigned int j=0; j<gen; j++) { //n° of generations
		unsigned int i = 0;
		while (i < N){
			if (i != N-1){
			x = rnd.Rannyu(0,1);
			if(x<0.4){
				n = rnd.Discrete(0,3);
				if (n==0)pop.PairPermutation(i);
				if (n==1)pop.Shift(i);
				if (n==2)pop.Permutation(i);
				if (n==3)pop.Inversion(i);
				i += 1;
			}
			else{
				pop.Crossover(i,i+1);
				i += 2;
			}
			}
			else{
				x = rnd.Rannyu(0,1);
				if(x<0.4){
					n = rnd.Discrete(0,3);
					if (n==0)pop.PairPermutation(i);
					if (n==1)pop.Shift(i);
					if (n==2)pop.Permutation(i);
					if (n==3)pop.Inversion(i);
					i+=1;
				}
			}
			//////////////
			/////MIGRATION

			if ( (j%16==0) and (i<N-3) ){
				for (unsigned int m = 0; m< 3; m++){
					for (unsigned int h = 0; h < cs; h++){
						migr_send[m*cs+h] = pop.GetPos(pop.GetMember(N-m-1),h);
					}
				}

				for (int k = 0; k < size; k++){
					if (rank == k) MPI_Isend(migr_send, 3 * cs, MPI_INTEGER, (k + shift) % size, k, MPI_COMM_WORLD, &req);	
					if (rank == (k + shift) % size)MPI_Recv(migr_rec, 3 * cs, MPI_INTEGER, k, k, MPI_COMM_WORLD, &stat);
				}

				for (unsigned int m = 0; m< 3; m++){
					// cout << "osl" << endl;
					for (unsigned int h = 0; h < cs; h++){
						// cout << "l" << endl;
						appo2.push_back(migr_rec[m*cs+h]);
						pop.SetSon(appo2,i+m);
					}
				}
			}			
		}
		if (j%16 == 0 ) cout << "n = " << j << endl;
		pop.SetPopulation(pop.GetSon());
		pop.Sorter();

		double accu = 0;
		for(unsigned int p = N/2; p < N; p++)accu += pop.NormL(pop.GetMember(p), 0);
		if (print.is_open()) print << 2.*accu/N << " " << pop.NormL(pop.GetMember(N-1), 0) << " " ;

	}
	delete[] migr_rec;
	delete[] migr_send;
	print.close();

/*
	if ( i % 10 == 0 )
		{
			//shift = 2*rnd.RannyuDiscr()-1;
			//if (rank == 0) shift = rnd.RannyuDiscr(1,4);
			//MPI_Bcast(&shift, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
			shift = 1;
			nmigr = nroutes*fractomigr;
			migrants1D_rec = new int[nmigr*ncities];
			for (int k = 0; k < nmigr; k++)
			{
				migrants.push_back(newpop[k]);
				newpop.erase(newpop.begin()+k);
			}
			migrants1D_send = resize<int>(migrants);
			//cout << "sendcheck" << rank << endl;
			for (int k = 0; k < size; k++)
			{
				if (rank == k)
					MPI_Isend(migrants1D_send, nmigr * ncities, MPI_INTEGER, (k + shift) % size, k, MPI_COMM_WORLD, &req);
				if (rank == (k + shift) % size)
					MPI_Recv(migrants1D_rec, nmigr * ncities, MPI_INTEGER, k, k, MPI_COMM_WORLD, &stat);
			}
			//cout << "postcheck" << rank << endl;			
			//if (rank == 1 ){double a = rnd.Rannyu();while (a > 0.0000000001) a = rnd.Rannyu();}
			migrants.clear();
			migrants = resize<int>(migrants1D_rec, nmigr * ncities, ncities);
			newpop.insert(newpop.begin(), migrants.begin(), migrants.end());
			//FittingSort(newpop, map);
			//newpop.assign(newpop.begin(), newpop.begin() + nroutes);
			delete[] migrants1D_rec;
			//cout << "check" << rank << endl;
			//cout << "doublecheck" << rank << endl;
			delete[] migrants1D_send;
			migrants.clear();
			//cout << "migrcheck" << rank << endl;
		}



	MPI_Gather(pop1D, N * cs, MPI_INTEGER, world1D, N * cs, MPI_INTEGER, 0, MPI_COMM_WORLD);
*/

/*
	for (unsigned int i=0; i<N; i++) {
		cout << "L_1 = " << pop.NormL(pop.GetMember(i), 0);
		cout << " ; L_2 = " << pop.NormL(pop.GetMember(i), 1) << endl;
	}
	*/
	/*
	unsigned int *pop1D = new unsigned int[N*cs];
	for (unsigned int q = 0; q<N;q++){
		for(unsigned int r = 0; r < cs; r ++){
			pop1D[cs*q+r] = pop.GetPos(pop.GetMember(q),r);
		}
	}
	int *world1D = new int[N*cs*size];
	cout << "?" << endl;
	MPI_Gather(pop1D, N * cs, MPI_INTEGER, world1D, N * cs, MPI_INTEGER, 0, MPI_COMM_WORLD);
	cout << "!" << endl;
	delete[] pop1D;
	*/

	ofstream write;
	write.open("citysquare" + to_string(rank) +".out");
	cout << " rank n° : " << rank << endl;
	for (unsigned int k = 0; k < cs; k++){
		//for (unsigned int i=0; i<Sim.GetN(); i++) 
		unsigned int i = pop.GetPos(pop.GetMember(N-1),k);
		cout << "i = " << i;
		cout << " , città n° " << k+1 << ": ( " << Sim.GetPosX(i-1) << " ; "<< Sim.GetPosY(i-1) <<  " )" << endl;
		if (write.is_open())write << Sim.GetPosX(i-1) << " " << Sim.GetPosY(i-1) << " " ;
		else cerr << "PROBLEM: Unable to open citysquare" << rank << ".out" << endl;
	}
	cout << "Total length = " << pop.NormL(pop.GetMember(N-1),0) << endl;
	write.close();

	/*
	unsigned int *migr_send_ = new unsigned int[cs * 3];
	unsigned int *migr_rec_ = new unsigned int[cs * 3];

	for (unsigned int j=0; j<5; j++) { //n° of generations
		unsigned int i = 0;
		while (i < N){
			if (i != N-1){
			x = rnd.Rannyu(0,1);
			if(x<0.4){
				n = rnd.Discrete(0,3);
				if (n==0)pop.PairPermutation(i);
				if (n==1)pop.Shift(i);
				if (n==2)pop.Permutation(i);
				if (n==3)pop.Inversion(i);
				i += 1;
			}
			else{
				pop.Crossover(i,i+1);
				i += 2;
			}
			}
			else{
				x = rnd.Rannyu(0,1);
				if(x<0.4){
					n = rnd.Discrete(0,3);
					if (n==0)pop.PairPermutation(i);
					if (n==1)pop.Shift(i);
					if (n==2)pop.Permutation(i);
					if (n==3)pop.Inversion(i);
					i+=1;
				}
			}
			//////////////
			/////MIGRATION

			if ( (j%16==0) and (i<N-3) ){
				for (unsigned int m = 0; m< 3; m++){
					for (unsigned int h = 0; h < cs; h++){
						migr_send_[m*cs+h] = pop.GetPos(pop.GetMember(N-m-1),h);
					}
				}

				for (int k = 0; k < size; k++){
					if (rank == k) MPI_Isend(migr_send_, 3 * cs, MPI_INTEGER, (k + shift) % size, k, MPI_COMM_WORLD, &req);	
					if (rank == (k + shift) % size)MPI_Recv(migr_rec_, 3 * cs, MPI_INTEGER, k, k, MPI_COMM_WORLD, &stat);
				}

				for (unsigned int m = 0; m< 3; m++){
					// cout << "osl" << endl;
					for (unsigned int h = 0; h < cs; h++){
						// cout << "l" << endl;
						appo2.push_back(migr_rec_[m*cs+h]);
						pop.SetSon(appo2,i+m);
					}
				}
			}			
		}
		if (j%1 == 0 ) cout << "*n = " << j << endl;
		pop.SetPopulation(pop.GetSon());
		pop.Sorter();

		double accu = 0;
		for(unsigned int p = N/2; p < N; p++)accu += pop.NormL(pop.GetMember(p), 0);
		if (print.is_open()) print << 2.*accu/N << " " << pop.NormL(pop.GetMember(N-1), 0) << " " ;

	}
	delete[] migr_rec_;
	delete[] migr_send_;
	*/
	print.close();
	
	MPI_Finalize();
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
