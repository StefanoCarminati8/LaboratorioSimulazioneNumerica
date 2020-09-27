#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"
#include "metodoblocchi.h"
#include "Integrale.h"
#include "funzioni.h"

using namespace std;

int main (int argc, char *argv[]){

  int nIntegrali = 10000;
  int nblk = 100;

  //Primo punto
  //Funzione da integrare
  Funzione* coseno = new Coseno();
  Integrale Integrale1(0,1, coseno);

  vector<double> integrali1(nIntegrali,0);
  for(int i=0; i<nIntegrali; i++){

    integrali1.at(i) = Integrale1.Media(nIntegrali);

  }

  Calcolo(nblk, integrali1, "risultati02.1.1.dat");

  //Secondo punto
  //Funzione da integrare
  Funzione* g = new gfunz();
  Integrale Integrale2(0,1, g);

  vector<double> integrali2(nIntegrali,0);
  for(int i=0; i<nIntegrali; i++){

    integrali2.at(i) = Integrale2.MediaRetta(nIntegrali);

  }

  Calcolo(nblk, integrali2, "risultati02.1.2.dat");

  delete coseno;
  delete g;

  return 0;
}
