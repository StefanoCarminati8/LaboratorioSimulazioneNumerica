#include <cmath>
#include "funzioni.h"

using namespace std;

//Implementazione delle funzioni da integrare

Funzione::~Funzione(){}

//Primo punto
double Coseno::Eval(double x) const {
	return (M_PI/2)*cos(x*(M_PI/2));
}
Coseno::~Coseno(){}

//Secondo punto
double gfunz::Eval(double x) const {
	return ((M_PI/2)*cos(x*(M_PI/2)))/(2*(1-x));
}

gfunz::~gfunz(){}
