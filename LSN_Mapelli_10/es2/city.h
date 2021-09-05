/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __city__
#define __city__
#include <vector>
#include <iostream>
#include "random.h"
using namespace std;

class city : public Random{

private:
  vector<vector<double>> position;
  unsigned int m_N;
  Random rand;
protected:

public:
  // constructors
  city();
  city(unsigned int, Random);
  // destructor
  ~city();
  // methods
  vector<double> GetPos(unsigned int)const;
  double GetPosX(unsigned int)const;
  double GetPosY(unsigned int)const;
  void SetPosX(unsigned int, double);
  void SetPosY(unsigned int, double);
  unsigned int GetN()const;
  void circle1D();
  void square2D();
  void print();

};

#endif // __city__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
