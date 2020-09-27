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
  int _stepfatti;
  int _stepbuoni;
  int _stepscartati;

public:

  Metropolis(Function *DistributionToBeSampled, vector <double> InitialPositions, double StepWidth, int nsteps, int dim);
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
  vector<double> raggi();

};

#endif
