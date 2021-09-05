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
#include <algorithm>
#include <vector>
#include "city.h"
#include "random.h"
#include "GAs.h"

using namespace std;

GAs :: GAs(){}

GAs :: GAs(unsigned int C, unsigned int R, city sim, Random rnd){
   m_C = C;
   m_R = R;
   Sim = sim;
   rand = rnd;

   for (unsigned int i=0; i<m_C; i++)member.push_back(0);
   for (unsigned int i=0; i<m_C; i++)member[i]=i+1;
   
   for (unsigned int i=0; i<m_R; i++){ 
      vector<unsigned int> zero (m_C);
      population.push_back(zero);   
   }
   
   population[0] = member;

   for (unsigned int i=1; i<m_R; i++){ 
      unsigned int j = 0;
      population[i][0] = 1;
      for (unsigned int n=2; n<m_C+1; n++){
         do{
            j = rand.Discrete(1,m_C-1); //pick a random place in the vector
            //cout << j << "-" ;
         }while (population[i][j] != 0);
         //cout << endl;
         population[i][j] = n;
      }
   }
   son = population;
   return;     
}

GAs :: ~GAs(){}

unsigned int GAs :: GetPos(vector<unsigned int> member, unsigned int N)const{
   return member[N];
}

vector<unsigned int> GAs ::GetMember (unsigned int n){
   return population[n];
}

vector<vector<unsigned int>> GAs ::GetPopulation (){
   return population;
}

vector<vector<unsigned int>> GAs ::GetSon (){
   return son;
}

void GAs :: SetSon (vector<unsigned int> member, unsigned int l){
   son[l] = member;
}

unsigned int GAs :: GetC()const{
   return m_C;
}
unsigned int GAs :: GetR()const{
   return m_R;
}

double GAs :: NormL( vector<unsigned int> membro, bool bit){
   //vector<unsigned int>::iterator it;
   double accu=0;
   double accu2=0;
   member = membro;
   for (unsigned int j=0; j<m_C; j++){
      unsigned int ci = member[j]-1;
      //cout << "member [" << j << "] = " << member[j] << " ; ci = " << ci << endl;
      //cout << "? = " << Sim.GetN() << "; " << m_C <<  endl;
      unsigned int cf = member[(j+1)%m_C]-1;
      //cout << "member [" << j+1 << "] = " << cf << endl;
      //cout << ci << "," << cf << endl;
      //cout << Sim.GetPosX(cf) << endl;
      accu += sqrt(  pow(Sim.GetPosX(cf)-Sim.GetPosX(ci), 2) + pow(Sim.GetPosY(cf)-Sim.GetPosY(ci), 2)   ) ;
      //cout << endl << "accu = " << accu << endl;
      //cout << j << "/" << m_C << endl;
      accu2 += pow(Sim.GetPosX(cf)-Sim.GetPosX(ci), 2) + pow(Sim.GetPosY(cf)-Sim.GetPosY(ci), 2) ;
      //cout << endl << "accu2 = " << accu2 << endl;
      //cout << j << "/" << m_C << endl;
   }
   if (bit == 0)return accu;
   if (bit == 1) return accu2;
   else{return 0;}
}

bool GAs :: Check (){
   int accu = 0;
   for (unsigned int i=0; i<m_R; i++){
      if (population[i][0] == 1){
         for (unsigned int j=1; j<m_C; j++){
            for (unsigned int k=j+1; k<m_C; k++){
               if (population[i][k] == population[i][j]){
                  cout << "error : repeating city ... " << population[i][k] << endl;
                  accu += 1;
                  return 1;
               }
               if (population[i][k] == 0){
                  cout << "error : city n°0 ... in chromosome " << i << " , gene " << k << endl;
                  accu += 1;
                  return 1;
               }

               else {accu +=0;}                  
            }
         }
      }
      else {
         cout << "error : the first city is not the n°1 ! " << endl;
         accu += 1;
         return 1;
      }
   }
   return accu;
}





/*
double Norm( vector<unsigned int> membro, bool bit, city Sim){
   double accu=0;
   double accu2=0;
   for (unsigned int j=0; j<Sim.GetN(); j++){
      unsigned int ci = membro[j]-1;
      unsigned int cf = membro[(j+1)%m_C]-1;
      accu += sqrt(pow(Sim.GetPosX(cf)-Sim.GetPosX(ci), 2) + pow(Sim.GetPosY(cf)-Sim.GetPosY(ci), 2)) ;
      accu2 += pow(Sim.GetPosX(cf)-Sim.GetPosX(ci), 2) + pow(Sim.GetPosY(cf)-Sim.GetPosY(ci), 2) ;
   }
   if (bit == 0)return accu;
   if (bit == 1) return accu2;
   else{return 0;}
}
}
*/

void GAs :: Sorter(){

   vector <double> l;
   for (unsigned int i = 0; i < m_R ; i++)l.push_back(NormL(GetMember(i),0));
   sort(l.begin(), l.end());
   vector <double> L = l;
   for (unsigned int i = 0; i < m_R ; i++)L[i] = l[m_R-1-i];

   vector <vector<unsigned int>> appo = population; 
   for (unsigned int i = 0; i < m_R ; i++){
      for (unsigned int j = 0; j < m_R ; j++){
         //cout << "norma [" << j << "] =" << NormL(population[j],1) << endl;
         //cout << "L[" << i<< "] = " << L[i] << endl;
         if (NormL(population[j],0) == L[i]){
            //cout << NormL(population[j],1) << endl;
            appo[i] = population[j];	
		      //for	(unsigned int j=0; j<m_C; j++)cout << appo[i][j] << " - " ;		            
		   }
      }
   }
    /*          
   for (unsigned int i=0; i<m_R; i++) {
		for	(unsigned int j=0; j<m_C; j++){
			cout << appo[i][j] << " - " ;		
		}
      cout << endl;
   }*/
   population = appo;
}
//sort (population.begin(), population.end(), Comparator);

/*
class mycomparatorfunctor {
public :
   bool operator () (vector <unsigned int> m1, vector <unsigned int> m2 ){return NormL(m1,0) < NormL(m2,0); ;};
} funct;
//bool mycomparatorfunction (vector <unsigned int> m1, vector <unsigned int> m2) {return NormL(m1,0) < NormL(m2,0) ; } ;
//sort (population[0], population[Sim.GetN()], Comparator(i,j)) ;
for (unsigned int i = 0; i < Sim.GetN(); i++){
   NormL(population[i],0);
}
//bool mycomparatorfunction (vector <unsigned int> m1, vector <unsigned int> m2) {return GAs.NormL(m1,0) < NormL(m2,0) ; } ;
sort( population[0], population[Sim.GetN()-1], Comparator() );
//sort (population[0], population[Sim.GetN()], [Comparator] (vector <unsigned int> i, vector <unsigned int> j) {Comparator(i,j) ; } ) ;
*/

unsigned int GAs :: Selection(){ //fitness selection
   double accu = 0;
   double accu2=0;
   double sum = 0;
   unsigned int answer = 0;
   vector <double> p;
   vector <double> p_sum;

   for (unsigned int i = 0; i< m_R ; i++)accu += NormL(population[i],1);
   
   for (unsigned int i = 0; i< m_R ; i++)p.push_back(1-(NormL(population[i],1)/accu));
   
   for (unsigned int i = 0; i< m_R ; i++)accu2 += p[i];

   for (unsigned int i = 0; i< m_R ; i++)p[i] /= accu2;

   for (unsigned int i = 0; i< m_R ; i++){
      sum += p[i];
      p_sum.push_back(sum);
   }
   double x = rand.Rannyu(0,1);
   //cout << endl << " x = " << x << endl;
   for (unsigned int i = 0; i< m_R ; i++){
      //cout << "i = " << i << endl;
      //cout << "p[i] = " << p_sum[i] << " , " << p_sum[i+1]<< endl;
      if( (x > p_sum[i]) & (x < p_sum[i+1]))answer=i;
      /*else {
         cout << "PROBLEM!" << endl;
         return 0;}*/
   }
   //cout << "selected n° : " << answer << endl;
   return answer;
}
unsigned int GAs :: SelectionFit (){ //fitness selection
   //Sorter();
   double fitness = 0;
   double sum = 0;
   double accu = 0;
   unsigned int answer = 0;
   vector <double> p;
   vector <double> p_sum;

   for (unsigned int i = 0; i< m_R ; i++){
      fitness = NormL(population[0],0) - NormL(population[i],0);
      p.push_back(fitness);
   }
   for (unsigned int i = 0; i< m_R ; i++)accu += p[i];
   for (unsigned int i = 0; i< m_R ; i++){
      sum += p[i];
      p_sum.push_back(sum/accu);
   }
    double x = rand.Rannyu(0,1);
   for (unsigned int i = 0; i< m_R ; i++){
      if( (x > p_sum[i]) & (x < p_sum[i+1]))answer=i;
   }
   return answer;
}

void GAs :: SetPopulation(vector<vector<unsigned int>> member){
   population = member;
}

/*bool GAs :: Twins(){
   int answer = 0;
   for (unsigned int i=0; i<m_R; i++){
      for (unsigned int j=i; j<m_R; j++){
               for (unsigned int k=0; i<m_C; i++)if (son[j][k]==son[i][k])answer++;
               //if (son[j] == son[i]) answer++;
         else return 0;
      }
   }

  return answer;
}*/




bool GAs :: Comparator (vector <unsigned int> m1, vector <unsigned int> m2){
   return NormL(m1,1) > NormL(m2,1);
}

void GAs :: pop (){
      int u = son.size();
      for (int k=0; k <u; k++){
         son.pop_back();
      }


   return;


}


///////////////////////
//GENETICS ALGORITHMS//
///////////////////////

unsigned int pbc (unsigned int, unsigned int); // periodic boundary condition without first position
unsigned int pbc (unsigned int i, unsigned int m_C){
   return ((i-1)%(m_C-1))+1;
   }

void GAs :: PairPermutation (/*vector <unsigned int> m1,*/unsigned int l){
   unsigned int n = SelectionFit();
   unsigned int i,j;
   member = population[n];
      i = rand.Discrete(1,m_C-1);
      do {
         j = rand.Discrete(1,m_C-1);
      } while (j == i);
   unsigned int appo = member[i];
   member[i] = member[j];
   member[j] = appo;
   //cout << "i = " << i << ", j = " << j << endl; 
   Check();
   son[l]=member;
}

void GAs :: Shift(/*vector <unsigned int> m1,*/ unsigned int l){ 

   unsigned int n = SelectionFit();
   unsigned int i,k,m; //??
   member = population[n];

   //unsigned int l=1;
   i = rand.Discrete(1,m_C-1); //first gene
   m = rand.Discrete(1,m_C-1); //how many genes
   k = rand.Discrete(1,m_C-2); //displacement
   //cout << "n = " << n << " ; i = " << i << " ; k = " << k << " ; m = " << m << endl;

   //cout << "n = " << n << " ; i = " << i << endl;
   vector <unsigned int> appo = member;
   //cout << " - = -";
   //for (unsigned int w = 0; w < m_C; w++) cout<< " * " << population[n][w] << " , ";

   for (unsigned int z = 0; z < m; z++){
      appo[pbc(i+k+z,m_C)] = member[pbc(i+z,m_C)];
      //cout <<" ... " << pbc(i+k+z,m_C) << " e " << pbc(i+z,m_C);
      //for (unsigned int w = 0; w < m_C; w++) cout << " * " << appo[w] << " , ";
      
      if(i+k+z-m_C+1!=i){
      //if ((pbc(i+z,m_C)<=pbc(i+k+z,m_C))&(pbc(i+z,m_C)>=)){
         //cout << "IF" << pbc(i+k+z,m_C);
         for (unsigned int w = 0; w < k; w++){
            //if(   (  pbc(i+w,m_C) <= pbc(i+k+z,m_C)) & (pbc(i+w,m_C)>= pbc(i+z,m_C)) )
           // {
               appo[pbc(i+w,m_C)] = member[pbc(i+w+z+1,m_C)] ;
               //cout <<" ,,, " << pbc(i+w,m_C) << " e " << pbc(i+w+z+1,m_C);

           }
      }
           else{
               break;
            }
         
      //}

     // else{
         
         //cout << "ELSE" << pbc(i+k+z,m_C);
         /*for (unsigned int w = 0; w < k; w++){
            appo[pbc(i+w+l,m_C)] = population[n][pbc(i+w+z+1,m_C)] ;
            cout <<" ,,, " << pbc(i+w+l,m_C) << " e " << pbc(i+w+z+1,m_C);
            
         }
         l++;
*/
         //break;
         /*
         for (unsigned int w = 0; w < k; w++){
            appo[pbc(i+w,m_C)] = population[n][pbc(i+w+z,m_C)] ;
            cout <<" ,,, " << pbc(i+w,m_C) << " e " << pbc(i+w+z,m_C);
         }*/
     // }
      
      
      //cout << endl << "z = " << z << " " ;
      //for (unsigned int w = 0; w < m_C; w++) cout << " * " << appo[w] << " , ";
      //appo[((i+k+j-1)%(m_C-1))+1] = population[n][((i-1+j)%(m_C-1))+1]; //j=0 : appo[8] = pop[9]; j = 1 : appo[9] = pop [1]
      //appo[((i+n+j-2)%(m_C-1))+1] = population[k][((i+j)%(m_C-1))+1]; //j=0 : appo[9] = pop[1]; j = 1 : appo[1] = pop [2]
   }

   member = appo;
   Check();
   //son.push_back(member);
   son[l]=member;

}
void GAs :: Permutation (/*vector <unsigned int> m1,*/ unsigned int l){
   
   unsigned int n = SelectionFit();
   member = population[n];
   unsigned int i,j,m; 
   m = rand.Discrete(2,int(floor((m_C-1)/2)-1));
   i = rand.Discrete(1,m_C-1);
   //cout << "m = " << m << " ; i = " << i << endl;
   do {
      j = rand.Discrete(1,m_C-1);
   } while (   min(  abs(int(i)-int(j))  ,    (((min(int(i),int(j))+int(m_C))-max(int(i),int(j)))%(int(m_C))) ) <=  int(m)  )  ;
   //cout << abs(int(i)-int(j)) << " " << (int(i)+int(j))%(int(m_C)-1)  << endl;

   //cout << "m = " << m << " ; i = " << i << " ; j = " << j  << endl;

   vector <unsigned int> appo = member;
   for (unsigned int k = 0; k < m; k++){
      appo[pbc(i+k,m_C)] = member[pbc(j+k,m_C)]; 
      appo[pbc(j+k,m_C)] = member[pbc(i+k,m_C)];
      //for (unsigned int w = 0; w < m_C; w++) cout << " * " << appo[w] << " , ";
      //cout<< endl;

      }
   member = appo;
   Check();
   //son.push_back(member);
   son[l]=member;


}

void GAs :: Inversion (/*vector <unsigned int> m1,*/unsigned int l){

   unsigned int n = SelectionFit();
   member = population[n];
   unsigned int i,m;
   i = rand.Discrete(1,m_C-1);
   m = rand.Discrete(2,m_C-1);
   //cout << "m = " << m << " ; i = " << i << endl;

   vector <unsigned int> appo = member;
   for (unsigned int k = 0; k < m; k++){
      appo[pbc(i+k,m_C)] = member[pbc(i+m-k-1,m_C)]; 
      //appo[pbc(i+m+k,m_C)] = population[n][pbc(i+k,m_C)]; 
      }
      //for (unsigned int w = 0; w < m_C; w++) cout << " * " << appo[w] << " , ";
     // cout<< endl;

   member = appo;
   Check();
   //son.push_back(member);
   son[l]=member;


}

void GAs :: Crossover (/*vector <unsigned int> m1,*/unsigned int l, unsigned int m){
      //k è il primo elemento del crossingover
   vector<unsigned int> appo1, appo2;
   unsigned int i,j,k;
   i = SelectionFit();
   member = population[i];
   j = SelectionFit();
   member1 = population[j];

   k = rand.Discrete(1,m_C-1);

   //cout << "k = " << k << " ; i = " << i << " ; j = " << j  << endl;

   for (unsigned int z=0; z<m_C-k; z++){
      //cout << ": " << population[i][k+z] << endl;
      appo1.push_back(member[k+z]);
      appo2.push_back(member1[k+z]);
   }

   for (unsigned int z=0; z<appo1.size(); z++){
      //cout << appo1[z] << " - " << appo2[z] << endl;
      //cout << "size " << appo1.size() << " " << appo2.size() << endl;
   }      
   //sort (appo1.begin(), appo1.end());
   //sort (appo2.begin(), appo2.end());

   for (unsigned int z=0; z<m_C-k; z++){
      //cout << Position(appo2, z) << endl;
      member[k+Position(appo2, z)] = appo1[Position(appo1,z)];
      //population[j][k+Position(appo1, z)] = population[j][appo2[Position(appo2,z)]]; 
      member1[k+Position(appo1, z)] = appo2[Position(appo2,z)];
   }
   Check();
   son[l]=member;
   son[m]=member1;

}

void GAs :: Nothing (/*vector <unsigned int> m1,unsigned int n*/){
   unsigned int n = SelectionFit();
   member = population[n];
   son.push_back(member);
}

unsigned int GAs :: Position (vector<unsigned int> member, unsigned int i){
   vector<unsigned int> appo = member;
   unsigned int min=0;
   unsigned int j=0;

   sort(appo.begin(), appo.end());
   min = appo[i];
   for (unsigned int k=0; k < member.size(); k++){
      if (member[k] == min){
         j=k;
      }
   }
   return j;
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
