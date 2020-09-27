#include "integrale.h"
#include "random.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

//Costruttore
Integrale::Integrale(double a, double b, Funzione * Funzione){

	_integranda = Funzione;
	//ordino gli estremi in ordine crescente e se non dovessero già essere ordinati ne tengo conto
	_a = min(a,b);
	_b = max(a,b);
	if ( _a>= _b){
		_Segno=-1;
	}
	else _Segno=1;

	Start();

}

//Distruttore
Integrale::~Integrale() {
	_rnd.SaveSeed();
}

void Integrale::Start(){

	// preparo il generatore di numeri casuali
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
	  	 _rnd.SetRandom(seed,p1,p2);
	 	}
	     }
	    input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

}

//Calcolo dell'integrale montecarlo usando numeri distribuiti uniformemente tra [_a,_b)
double Integrale::Media(int N){

	_sommaparziale = 0;
	double media = 0;
	double intervallo = _b - _a;

	for (int i=1; i <= N; i++){
		double x = _rnd.Rannyu(_a,_b);
		double f = _integranda -> Eval(x);
		_sommaparziale += f;
	}

	media = (_sommaparziale/double(N));
	Integrale = _Segno * media * intervallo;

	return Integrale;

};

double Integrale::MediaRetta(int N){

	_sommaparziale = 0;
	double media = 0;
	double intervallo = _b - _a;;

	for (int i=1; i <= N; i++){
		double x = _rnd.Retta();
		double f = _integranda -> Eval(x);
		_sommaparziale += f;
	}

	media = (_sommaparziale/double(N));
	Integrale = _Segno * media * intervallo;

	return Integrale;

};

//restituisce il valore dell'integrale
double Integrale::GetResult() const{
	return Integrale;
};
