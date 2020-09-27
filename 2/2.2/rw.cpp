#include "rw.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"
#include <string>


using namespace std;

//Costruttore
RW::RW(double LatticeStep){

	Start();
	//Inizializza la posizione nell'origine
	_x = 0;
	_y = 0;
	_z = 0;

	_LatticeStep = LatticeStep;
}

//Distruttore
RW ::~RW(){
	_rnd.SaveSeed();
}

void RW::Start(){

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

//Metodo per restartare il Random Walk reimpostando tutto a zero
void RW::Restart(){
	_x = 0;
	_y = 0;
	_z = 0;

}

//Metodo che ggiorna le posizioni e le distanze dopo un singolo passo in un reticolo cubico tridimensionale
void RW::CubicLatticeStep(){

	//Estraggo un intero casuale tra 0 e 5 e assegno ad ognuno una direzione arbitraria
	int p = int(_rnd.Rannyu(0,6));

	//Aggiorno la posizione
	if(p==0) _x++;
	if(p==1) _y++;
	if(p==2) _z++;
	if(p==3) _x--;
	if(p==4) _y--;
	if(p==5) _z--;

}

//Metodo che mi permette di effettuare NumberOfStep passi nel reticolo cubico
double RW::DiffusionLattice(int NumberOfSteps){

	for(int i=0; i<NumberOfSteps; i++) CubicLatticeStep();

	double d = _LatticeStep * (sqrt((pow(_x,2)+pow(_y,2)+pow(_z,2))));
	Restart();

	return d;
}

//Metodo che aggiorna le posizioni e le distanze dopo un singolo passo in direzione casuale nello spazio tridimensionale
void RW::IsotropicStep(){

	//Campiono uniformemente l'angolo solido
	double theta = _rnd.Rannyu(0,M_PI);
	double phi = _rnd.Rannyu(0,2*M_PI);

	//Aggiorno la posizione
	_x += _LatticeStep * sin(theta)*cos(phi);
	_y += _LatticeStep * sin(theta)*sin(phi);
	_z += _LatticeStep * cos(phi);

}

//Metodo che mi permette di effettuare NumberOfSteps passi nello spazio tridimensionale in una direzione casuale
double RW::IsotropicDiffusion(int NumberOfSteps){

	for(int i=0; i<NumberOfSteps; i++) IsotropicStep();

	double d= _LatticeStep * (sqrt((pow(_x,2)+pow(_y,2)+pow(_z,2))));
	Restart();

	return d;
	}
