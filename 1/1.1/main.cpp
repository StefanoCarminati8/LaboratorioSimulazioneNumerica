#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "metodoblocchi.h"

using namespace std;

int main() {

  Random rnd;
  SetRandomGenerator(&rnd);

  int N=100000; //numero dati
  int L=1000;  //numero blocchi
  int F=100;  //divisione dell'intervallo [0,1) es.3
  int M=1000000;  //numeri generati es.3 (10^4 per intervallo)

  vector<double> dati1(N,0);
  vector<double> dati2(N,0);
  vector<double> dati3(M,0);
  vector<double> chiquadro(F,0);

  for(int i=0; i<N; i++){
    dati1.at(i) = rnd.Rannyu();
  }

  for(int j=0; j<N; j++){
    dati2.at(j) = pow((rnd.Rannyu() - 0.5), 2);
  }

  for(int k=0; k<M; k++){
    dati3.at(k) = rnd.Rannyu();
  }

  Calcolo(L, dati1, "ris.01.1.1.dat");
  Calcolo(L, dati2, "ris.01.1.2.dat");
  ChiQuadro(F,M, dati3, "ris.01.1.3.dat");

  rnd.SaveSeed();
  return 0;
}
