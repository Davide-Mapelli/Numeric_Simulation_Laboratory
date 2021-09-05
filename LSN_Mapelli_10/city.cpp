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
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <algorithm>
#include "random.h"
#include "city.h"

using namespace std;

city :: city(){}

city :: city(unsigned int N, Random rnd){
   m_N = N;
   rand = rnd;
   //position  = new (double [m_N])[2];
   for (unsigned int i=0; i<m_N; i++){ 
      vector<double> zero = {0,0};
      position.push_back(zero);
   }
}

city :: ~city(){
//delete [] position;
}
unsigned int city :: GetN()const{
   return m_N;
}

vector<double> city :: GetPos(unsigned int i)const{
   return position[i];
}

double city :: GetPosX(unsigned int i)const{
   return position[i][0];
}
double city :: GetPosY(unsigned int i)const{
   return position[i][1];
}

/////////////////////////////////////////////////////////////
void city :: circle1D(){
   //rand = rnd;
   for (unsigned int i=0; i<m_N; i++){
      double theta = rand.Rannyu(0,2*M_PI);
      vector<double> x_y = {cos(theta),sin(theta)};
      position.insert(position.begin()+i,x_y);
   }
   return;
}

void city :: square2D(){
   //rand = rnd;
   for (unsigned int i=0; i<m_N; i++){
      double x = rand.Rannyu(-1,+1);
      double y = rand.Rannyu(-1,+1);
      vector<double> xy = {x,y};
      position.insert(position.begin()+i,xy);
   }    
   return;
}

void city :: print(){
   ofstream write;
	write.open("city.out");
	if (write.is_open())for (unsigned int i=0; i< m_N; i++)write << GetPosX(i) << " " << GetPosY(i) << " " ;
	else cerr << "PROBLEM: Unable to open random.out" << endl;
	write.close();

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
