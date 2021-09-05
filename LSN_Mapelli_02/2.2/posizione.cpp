#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "posizione.h"


using namespace std;
  // costruttori
  Posizione::Posizione(){
	m_x=0.;
	m_y=0.;
	m_z=0.;
			}
  Posizione::Posizione(double x, double y, double z){
	m_x=x;
	m_y=y;
	m_z=z;
						}
  // distruttore
  Posizione::~Posizione() {}
  // metodi
  double Posizione::getX() const{return m_x;}       // Coordinate cartesiane
  double Posizione::getY() const{return m_y;} 
  double Posizione::getZ() const{return m_z;} 

  void Posizione::setX(double var) {m_x = var;}       // Coordinate cartesiane
  void Posizione::setY(double var) {m_y = var;} 
  void Posizione::setZ(double var) {m_z = var;} 

  double Posizione::getR() const{return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);}       // Coordinate sferiche
  double Posizione::getPhi() const{return atan2(m_y,m_x);} 
  double Posizione::getTheta() const{return acos(m_z/getR());} 
  double Posizione::getRho() const{return sqrt(m_x*m_x+m_y*m_y);}     // raggio delle coordinate cilindriche
  double Posizione::Distanza(const Posizione& b) const {	// distanza da un altro punto
	return sqrt(	pow(getX()-b.getX(),2)+pow(getY()-b.getY(),2)+pow(getZ()-b.getZ(),2)	);
					
}
