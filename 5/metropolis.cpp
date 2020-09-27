
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstring>
#include "random.h"
#include "funzioni.h"
#include "metodoblocchi.h"
#include "metropolis.h"

using namespace std;

//Costruttore
Metropolis::Metropolis(Function *DistributionToBeSampled, vector <double> InitialPositions, double StepWidth, int nsteps, int dim): _NextPositions(dim){

  SetRandomGenerator(&_rnd);

  _DistributionToBeSampled = DistributionToBeSampled;
  _StepWidth = StepWidth;
  _nsteps = nsteps;
  _dim = dim;
  _stepfatti = 0;
  _stepbuoni = 0;
  _stepscartati = 0;
  _Positions.push_back(InitialPositions);

}

//Distruttore
Metropolis::~Metropolis(){
  _rnd.SaveSeed();
}

//Metodo che riinizializza tutto a zero per riiniziare il metropolis dall'ultima posizione in memoria
void Metropolis::Restart(){

  _stepfatti = 0;
  _stepbuoni = 0;
  _stepscartati = 0;

  _Positions.erase(_Positions.begin(), _Positions.end() -1)  ;

}

//Metodo per equilibrare il sistema
void Metropolis::Equilibrate(int nstepeq, const char * trial){

	int temp = _nsteps;
	_nsteps = nstepeq;

	//A seconda di quale distribuzione di trial sto usando utilizzo Run diversi
	if(strcmp(trial,"uniform") == 0){
		RunUniform();
	}
	else if(strcmp(trial,"gauss") == 0){
		RunGauss();
	}

	Restart();
	_nsteps = temp;

}

//Possibile passo uniforme
void Metropolis::NextUniform(void){

  for(int i=0; i<_dim; i++){
    _NextPositions.at(i) = _rnd.Rannyu(-_StepWidth, _StepWidth) + _Positions[_stepfatti][i];
  }

}

//Possibile passo gaussiano
void Metropolis::NextGauss(void){

  for(int i=0; i<_dim; i++){
    _NextPositions.at(i) = _rnd.Gauss(0, _StepWidth) + _Positions[_stepfatti][i];
  }

}

//Probabilità di accettazione
bool Metropolis::Acceptance(){

	//Distribuzione di probabilità valutata nel punto in cui mi trovo
	double actualProbability = _DistributionToBeSampled-> Eval(_Positions[_stepfatti]);
	//Distribuzione di probabilità valutata nel possibile nuovo punto
	double newProbability = _DistributionToBeSampled-> Eval(_NextPositions);
	//Rapporto tra le distribuizioni di probabilità
	double probabilityRatio = newProbability/actualProbability;
	//Probabilità di accettazione
	double acceptance = min(1., probabilityRatio);

	if(acceptance >= 1){
		return true;
	}

	else{
		double a = _rnd.Rannyu();
		if( a <= acceptance){
			return true;
		}
		return false;
	}

}

void Metropolis::UniformStep(){

	NextUniform();
	bool acceptance = Acceptance();

	if(acceptance == true){

    _Positions.push_back(_NextPositions);
    _stepbuoni++;

	}
	else{
		vector<double> temp(_dim);	//Vettore di appoggio
    	_stepscartati++;
		for(int i=0; i<_dim; i++){
			temp.at(i) = _Positions[_stepfatti][i];
		}
		_Positions.push_back(temp);

	}

	//Incremento del numero totale di steps fatti
	_stepfatti ++;

}

void Metropolis::GaussStep(){

	NextGauss();
	bool acceptance = Acceptance();

	if(acceptance == true){

    _Positions.push_back(_NextPositions);
    _stepbuoni++;

	}
	else{
		vector<double> appoggio(_dim);
    	_stepscartati++;
		for(int i=0; i<_dim; i++){
			appoggio.at(i) = _Positions[_stepfatti][i];
		}
		_Positions.push_back(appoggio);

	}

	//Incremento del numero totale di steps fatti
	_stepfatti ++;

}

//UnifromStep ripetuto nsteps volte
void Metropolis::RunUniform(){

	for(int i=0; i < _nsteps; i++){
		UniformStep();
	}

}

//GaussStep ripetuto _stepbuoni volte
void Metropolis::RunGauss(){

	for(int i=0; i < _nsteps; i++){
		GaussStep();
	}

}

//Metodo per controllare la regola empirica del 50%
void Metropolis::AcceptanceRate(){
	cout << endl << "Accepted steps = " << _stepbuoni << endl;
  	cout << endl << "Refused steps = " << _stepscartati << endl;
	cout << endl << "ACCEPTANCE RATE: " << ((double)_stepbuoni/(double)_nsteps) * 100. << "%" << endl;
}

//Metodo per stampare su file i punti campionati
void Metropolis::Print(const char* filename){

	ofstream Results;
	Results.open(filename);

	for(int j = 0; j < _nsteps+1; j++){
		for(int i=0; i <_dim; i++){
			Results <<  _Positions[j][i] << " " ;
		}
		Results << endl;
	}

  Results.close();
}

//Metodo legato all'esercizio 05.1
//Metodo per ottenere i raggi dei punti (x,y,z) campionati
vector<double> Metropolis::raggi(){

	vector<double> raggi(_nsteps+1,0);

	for(int j = 0; j < _nsteps+1; j++){
		for(int i=0; i < _dim; i++){
			raggi.at(j) += pow( _Positions[j][i], 2);
		}
	raggi.at(j) = sqrt(raggi.at(j));
	}

return raggi;
}
