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
Metropolis::Metropolis(Function *DistributionToBeSampled, vector <double> posIniziali, double StepWidth, int nsteps, int dim, int MCSteps, int nbins): _NextPositions(dim), _Psi(nbins, vector<double>(MCSteps)){
	//Random _rnd;
	//SetRandomGenerator(&_rnd);

	_DistributionToBeSampled = DistributionToBeSampled;
	_StepWidth = StepWidth;
	_nsteps = nsteps;
	_dim = dim;
	_StepsDone = 0;
	_StepAccepted = 0;
	_StepRefused = 0;
	_Positions.push_back(posIniziali);
	_Extremes = 3;							//Dato strettamento relativo all'esericizio 08.1
	_nbins = nbins;
	_NumberOfIteractions = MCSteps;

	//Inizializzazione matrice per realizzazione degli istogrammi
	for(int j=0; j<_NumberOfIteractions; j++){
		for(int i=0; i<_nbins;i++){
			_Psi[i][j]=0;
		}
	}
}

//Distruttore
Metropolis::~Metropolis(){
  _rnd.SaveSeed();
}

//Metodo che riinizializza tutto a zero per riiniziare il metropolis dall'ultima posizione in memoria
void Metropolis::Restart(){

	_StepsDone = 0;
	_StepAccepted = 0;
	_StepRefused = 0;

	_Positions.erase(_Positions.begin(), _Positions.end() -1);

}

//Metodo per equilibrare il sistema
void Metropolis::Equilibrate(int NumberOfEquilibrationSteps, const char * trial){

	int temp = _nsteps;
	_nsteps = NumberOfEquilibrationSteps;

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
    _NextPositions.at(i) = _rnd.Rannyu(-_StepWidth, _StepWidth) + _Positions[_StepsDone][i];
  }

}

//Possibile passo gaussiano
void Metropolis::NextGauss(void){

  for(int i=0; i<_dim; i++){
    _NextPositions.at(i) = _rnd.Gauss(0, _StepWidth) + _Positions[_StepsDone][i];
  }

}

//Probabilità di accettazione
bool Metropolis::Acceptance(){

	//Distribuzione di probabilità valutata nel punto in cui mi trovo
	double actualProbability = _DistributionToBeSampled-> Eval(_Positions[_StepsDone]);
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
    _StepAccepted++;

	}
	else{
		vector<double> temp(_dim);	//Vettore di appoggio
    	_StepRefused++;
		for(int i=0; i<_dim; i++){
			temp.at(i) = _Positions[_StepsDone][i];
		}
		_Positions.push_back(temp);

	}

	//Incremento del numero totale di steps fatti
	_StepsDone ++;

}

void Metropolis::GaussStep(){

	NextGauss();
	bool acceptance = Acceptance();

	if(acceptance == true){

    _Positions.push_back(_NextPositions);
    _StepAccepted++;

	}
	else{
		vector<double> appoggio(_dim);
    	_StepRefused++;
		for(int i=0; i<_dim; i++){
			appoggio.at(i) = _Positions[_StepsDone][i];
		}
		_Positions.push_back(appoggio);

	}

	//Incremento del numero totale di steps fatti
	_StepsDone ++;

}

//UnifromStep ripetuto nsteps volte
void Metropolis::RunUniform(){

	for(int i=0; i < _nsteps; i++){
		UniformStep();
	}

}

//GaussStep ripetuto _StepAccepted volte
void Metropolis::RunGauss(){

	for(int i=0; i < _nsteps; i++){
		GaussStep();
	}

}

//Metodo per controllare la regola empirica del 50%
void Metropolis::AcceptanceRate(){
	cout << endl << "----------------------------------------" << endl;
	cout << "Accepted steps = " << _StepAccepted << endl;
  	cout << "Refused steps = " << _StepRefused << endl;
	cout << "ACCEPTANCE RATE: " << ((double)_StepAccepted/(double)_nsteps) * 100. << "%" << endl;
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
vector<double> Metropolis::Rays(){

	vector<double> rays(_nsteps+1,0);

	for(int j = 0; j < _nsteps+1; j++){
		for(int i=0; i < _dim; i++){
			rays.at(j) += pow( _Positions[j][i], 2);
		}
	rays.at(j) = sqrt(rays.at(j));
	}

return rays;
}

//Metodi legati all'esercizio 08.1
//Metodo della media
double Metropolis::Integrate(Function *hamiltonian){

	double sum = 0;

	for (int i=0; i < _nsteps; i++){
		double f = hamiltonian -> Eval(_Positions.at(i));
		sum += f;
	}
	return (sum/double(_nsteps));
}

void Metropolis::Histogram(int j){

	double dr = (2*_Extremes)/(double)_nbins;
	int index;
	double dato;

	for(int i=0; i<_StepsDone; i++){
		dato = _Positions[i][0];
		if( abs(dato) < _Extremes){
			index = int((_Extremes + dato) / dr);
			_Psi[index][j] += 1;
		}
	}
  	//Normalizzazione istogramma
  	double sum = 0;
  	for(int i = 0; i < _Psi.size(); i++){
		sum += _Psi[i][j];
	}

  	for(int i = 0; i < _Psi.size(); i++){
    	_Psi[i][j] = _Psi[i][j] / (dr * sum);
  	}
}

void Metropolis::AnalysisHistogram(int NumberOfBlocks, const char* OutputFile){

	double dr = (2*_Extremes)/(double)_nbins;
	AnalysisMatrix(NumberOfBlocks, _Psi , dr, _Extremes, OutputFile);

}
