#ifndef _rw_h
#define _rw_h
#include "random.h"

class RW{

private:

  Random _rnd;
  //Coordinate
  double _x;
  double _y;
  double _z;
  //Passo reticolare
  double _LatticeStep;

public:

  RW(double LatticeStep);
  ~RW();

  void Start();
  void Restart();
	void CubicLatticeStep();
	double DiffusionLattice(int NumberOfSteps);
  void IsotropicStep();
	double IsotropicDiffusion(int NumberOfSteps); //permette di effettuare N passi nello spazio tridimensionale in una direzione casuale
};

#endif
