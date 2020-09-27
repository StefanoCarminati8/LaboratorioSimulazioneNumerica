#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "metodoblocchi.h"

using namespace std;

int main() {

  Random rnd;
  SetRandomGenerator(&rnd);

  ofstream scrivi;
  scrivi.open("ris01.2.dat");

  int N[4]={1,2,10,100}; //dati per realizzazione
  int V=10000;  //realizzazioni
  double lambda=1.; //valori da utilizzare nell'esercizio
  double mu=0.;
  double gamma=1.;

  vector<int> numdati (N, N + sizeof(N) / sizeof(int)); //per poter usare i valori dei dati

  vector<double> medieunif(V,0);
  vector<double> mediesp(V,0);
  vector<double> medielor(V,0);

  for(int i=0; i<numdati.size(); i++){

    vector<double> uni(numdati.at(i));
    vector<double> esp(numdati.at(i));
    vector<double> lor(numdati.at(i));

    for(int j=0; j<V; j++){
      for(int s=0; s<numdati.at(i); s++){
        uni.at(s)=rnd.Rannyu();
        esp.at(s)=rnd.Exp(lambda);
        lor.at(s)=rnd.Lorentz(mu, gamma);
      }
      medieunif.at(j) = Media(uni);
      mediesp.at(j) = Media(esp);
      medielor.at(j) = Media(lor);

      if (scrivi.is_open()) scrivi << setprecision(12) << medieunif.at(j) << " " << setprecision(12) << mediesp.at(j) << " " << setprecision(12) << medielor.at(j) << endl;
      else cerr << "Errore nell'aprire file dati" << endl;
    }
  }
  scrivi.close();
  rnd.SaveSeed();
  return 0;

}
