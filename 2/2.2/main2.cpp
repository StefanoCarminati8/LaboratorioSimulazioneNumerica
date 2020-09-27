#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"
#include "rw.h"
#include "metodoblocchi.h"

using namespace std;

int main (int argc, char *argv[]){

  int N = 1000000;   //Numero di Random Walk totali
  int nstep = 100;            //Numero di passi del Random Walk e coincide con il numero di blocchi che utilizzo per il data blocking

  //RANDOM WALK DISCRETO
  RW randomwalk1(1);                  //Random walk di passo reticolare 1

  int RWsperstep = int(N/nstep);    //Numero di random walk che voglio simulare per ogni passo
  vector<double> d1(N,0);          //Vettore delle distanze al quarato

  for(int i=0; i<nstep; i++){
		for(int j=0; j<RWsperstep; j++){
      int n = j +i*RWsperstep; //indice distanza n-esima: per ogni num di passi avrò 10^4 random walks
			d1[n] = (double) pow(randomwalk1.DiffusionLattice(i) , 2);
		}
	}
  Calcolo(nstep, d1, "risultati02.2.1.dat");

  //RANDOM WALK CONTINUO
  RW randomwalk2(1);   //Random walk di passo 1

  vector<double> d2(N,0);

  for(int i=0; i<nstep; i++){
		for(int j=0; j<RWsperstep; j++){
      int n = j +i*RWsperstep;
			d2[n] = (double) pow(randomwalk2.IsotropicDiffusion(i) , 2);
		}
	}
  Calcolo(nstep, d2, "risultati02.2.2.dat");

  return 0;
}
