#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstring>
#include "random.h"
#include "metodoblocchi.h"
#include "funzioni.h"
#include "metropolis.h"

using namespace std;

int main (int argc, char *argv[]){

  int dim = 1;
  double metropolisStep = 2.6;
  int nsteps = 10000;
  int nintegrals = 10000;
  int nblk = 100;
  int nbins = 200;
  int neqstep= 100;

  WaveFunction* PsiTrial = new WaveFunction();
  //Impongo i parametri ottimali ottenuti da Variational MonteCarlo
  PsiTrial->SetMi(0.795);
  PsiTrial->SetSigma(0.61);
  Function* H = new Hamiltonian(PsiTrial);

	vector<double> InitialPositions(dim,0);
	InitialPositions.at(0) = 0;
  Metropolis walker(PsiTrial, InitialPositions, metropolisStep, nsteps, dim, nintegrals, nbins);
  vector<double> results_integration(nintegrals,0);

  walker.Equilibrate(neqstep, "uniform");
  for(int i=0; i<nintegrals; i++){

    walker.RunUniform();
    results_integration.at(i) = walker.Integrate(H);
    if(i%10000){
      walker.AcceptanceRate();
    }
    walker.Print("psi.out");
    walker.Histogram(i);
    walker.Restart();
  }

  walker.AnalysisHistogram(nblk,"configurations.out");                //DataBlocking applicato all'instogramma
  Calcolo(nblk, results_integration,"hamiltonian_results.out");      //DataBlocking degli integrali

  delete PsiTrial;
  delete H;
  return 0;

}
