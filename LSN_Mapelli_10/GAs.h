/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __GAs__
#define __GAs__

#include <vector>
#include <iostream>
#include "random.h"
#include "city.h"

class GAs : public city{

private:
vector<vector<unsigned int>> population;
vector<unsigned int> member;
vector<unsigned int> member1;
vector<unsigned int> path;
vector<unsigned int> newpath;
vector <vector<unsigned int>> son;

unsigned int m_C, m_R; //columns and rows
//Random rand;
city Sim;
Random rand;
protected:

public:
  // constructors
  GAs();
  GAs(unsigned int C, unsigned int R, city Sim, Random rand);
  // destructor
  ~GAs();
  // methods
  unsigned int GetPos(vector<unsigned int> member, unsigned int)const;
  vector<unsigned int> GetMember(unsigned int);
  vector<unsigned int> GetPath();
  vector<unsigned int> GetNewPath();

  void SetPath(vector<unsigned int> member);
  void SetNewPath(vector<unsigned int> member);

  vector<vector<unsigned int>> GetPopulation();
  vector<vector<unsigned int>> GetSon();
  unsigned int GetC()const;
  unsigned int GetR()const;

  double NormL(vector<unsigned int> member, bool bit);
  bool Check ();
  void Sorter ();  
  bool Comparator (vector <unsigned int> m1, vector <unsigned int> m2);
  unsigned int Selection();
  unsigned int SelectionFit();
  //void Results();
  unsigned int Position(vector<unsigned int> member, unsigned int z);
  void SetPopulation(vector<vector<unsigned int>> member);
  //bool Twins();
  void pop();

//GENETICS
  void PairPermutation ();
  void Shift();
  void Permutation();
  void Inversion();
  void Crossover(unsigned int i, unsigned int j);
  void Nothing();

};

#endif // __GAs__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
