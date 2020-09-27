#ifndef _metropolis_h
#define _metropolis_h

#include "random.h"
#include "funzioni.h"
#include <vector>

using namespace std;

class Metropolis{

private:

  Random _rnd;
  Function *_DistributionToBeSampled;
  int _dim;
  double _StepWidth;
  int _nsteps;
  vector <vector <double> > _Positions;
  vector <double> _NextPositions;
  int _StepsDone;
  int _StepAccepted;
  int _StepRefused;

  //Istogramma
  double _Extremes; //Estremi per il calcolo istogramma
	int _nbins;
	vector<vector<double>> _Psi; //Serve per l'istogramma delle psi
	int _NumberOfIteractions; //Numero di istogrammi su cui mediare

public:

  Metropolis(Function *DistributionToBeSampled, vector <double> posIniziali, double StepWidth, int nsteps, int dim, int MCSteps, int nbins);
  ~Metropolis();

  void Restart();
  void Equilibrate(int nstepeq, const char * trial);
  void RunUniform(void);
  void RunGauss(void);
  void NextUniform(void);    //Distribuzione trial == uniforme
  void NextGauss(void);      //Distribuzione trial == gaussiana
  bool Acceptance(void);     //Probalità di accettazione
  void UniformStep(void);
  void GaussStep(void);
  void AcceptanceRate(void);
  void Print(const char*);
  double Integrate(Function *hamiltonian);
  void Histogram(int index); //Aggiornamento dell'istogramma della funzione d'onda di trial
	void AnalysisHistogram(int nblk, const char* OutputFile);//Data blocking sull' istogramma di psi

  vector<double> Rays();
  vector<vector <double>> *GetSampleWafefunction(void);
};

#endif
