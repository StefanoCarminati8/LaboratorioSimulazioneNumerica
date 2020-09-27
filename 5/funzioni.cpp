#include <cmath>
#include <vector>
#include "funzioni.h"

using namespace std;

//Raggio di bohr
#define a0 1

// Le distribuzione da campionare
//Le distanze sono espresse in unità del raggio di Bohr!
Function::~Function(){}

double Uno_s::Eval(vector <double> pos) const {

	double r = sqrt(pow(pos.at(0),2) + pow(pos.at(1),2) + pow(pos.at(2),2));			//Trasformo in coordinate polari
	return pow(pow(a0 , - 3./2.) * exp(-r/a0)/sqrt(M_PI), 2);							//é come se fosse il modulo quadro tanto sono reali
}
Uno_s::~Uno_s(){}

double Due_p::Eval(vector <double> pos) const {

	double r = sqrt(pow(pos.at(0),2) + pow(pos.at(1),2) + pow(pos.at(2),2));			//Trasformo in coordinate polari
	return pow((pow(a0, - 5./2.) / 8.) * sqrt(2. / M_PI) * (exp(-r/(2.*a0)) * pos.at(2) ), 2);

}
Due_p::~Due_p(){}
