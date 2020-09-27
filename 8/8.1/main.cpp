#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "metodoblocchi.h"
#include "funzioni.h"
#include "metropolis.h"

using namespace std;

int main (int argc, char *argv[]){

  int dim = 1;
  double metrostep = 2.6;
  int steps = 1000;
  int nintegrali = 10000;
  int nstepeq = 1000;
  int nblk = 100;       //Inutile per questo punto
  int nbins = 200;         //Inutile per questo punto

  WaveFunction* PsiTrial = new WaveFunction();
  Function* H = new Hamiltonian(PsiTrial);

	vector<double> posIniziali(dim,0);
	posIniziali.at(0) = 0;

  Metropolis walker(PsiTrial, posIniziali, metrostep, steps, dim, nintegrali, nbins);
  vector<double> risultatiIntegrazione(nintegrali,0);
  ofstream parametri;
  parametri.open("optimization.parameter.prova.out");

  walker.Equilibrate(nstepeq, "uniform");
  for(double mi=0.795; mi<0.81; mi+=0.001){
    cout << "Mi: " <<mi << endl;
    PsiTrial->SetMi(mi);
    for(double sigma=0.61; sigma<0.62; sigma+=0.001){
      cout << "Sigma: " << sigma << endl;
      PsiTrial->SetSigma(sigma);
      for(int i=0; i<nintegrali; i++){
        walker.RunUniform();
        risultatiIntegrazione.at(i) = walker.Integrate(H);
        walker.Restart();
      }
      parametri << mi << "  " << sigma << "  " << Media(risultatiIntegrazione) << endl;
    }
  }
  parametri.close();

  delete PsiTrial;
  delete H;
  return 0;
}
