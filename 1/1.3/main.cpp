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

  double ago=0.; //posizione del centro dell'ago
  double d=7.; //distanza tra linee
  double L=3.; //lughezza ago
  double angolo=0.; //angolo a cui atterra l'ago
  double sup,inf=0.; //estremi superiore e inferiore dell'ago
  double ntiri=10000; //numero di volte che viene lanciato l'ago
  double nhit=0.; //numero di volte che l'ago attraversa una linea
  double prove=1000; //numero di volte che viene ripetuto l'esperimento

  vector<double> pi(prove,0);//vettore che conterrà i valori calcolati di pigreco

  for(int t=0; t<prove; t++){ //ciclo sul numero di esperimenti
    nhit=0.; //reset del numero dei risultati
    for(int i=0; i<ntiri; i++){
      angolo = rnd.Rannyu(0.,180.); //genero l'angolo casuale
      ago = rnd.Rannyu(-d*0.5,d*0.5);
      /*
      qui sopra, genero la posizione casuale dove cade il centro dell'ago
      compresa tra + e - la metà della distanza fra due linee. In questo
      modo, la linea si troverà al punto x=0 e sarà possibile valutare facilmente
      se l'ago si sovrapporà a una linea o meno cercando il valore dei suoi
      estremi e controllando che siano discordi.
      */
      inf = ago - sin(angolo*M_PI/180.)*L*0.5;
      sup = ago + sin(angolo*M_PI/180.)*L*0.5;

      if((inf <=0 && sup >=0) || (inf>=0 && sup<=0)){ //controllo se l'ago sia
                                                      //sopra una riga o no
        nhit = nhit+1;
      }
    }
    pi.at(t) = (2*L*ntiri)/(nhit*d); //calcolo il valore di pigreco ottenuto e
                                     //lo salvo in un vettore
  }

  Calcolo(100, pi, "ris.01.3.dat"); //utilizzo la funzione calcolo del metodoblocchi.cpp
                                    //per salvare in un file media e dev.std. dei
                                    //valori di pigreco ottenuti

  rnd.SaveSeed();
  return 0;
}
