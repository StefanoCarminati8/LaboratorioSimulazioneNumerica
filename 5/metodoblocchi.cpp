#include "metodoblocchi.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

//Somma di elementi di un vettore:
double Somma(vector<double> vettore, int inizio, int fine){
  double somma=0;
  for(int i=inizio; i<fine; i++){
    somma += vettore.at(i);
  }
  return somma;
}

//Media dei dati contenuti in un vettore:
double Media(vector<double> vettore){
  double somma=0;
  for (int i=0; i<vettore.size(); i++){
    somma += vettore.at(i);
  }
  return somma/vettore.size();
}

//Media dei quadrati dei dati contenuti in un vettore:
double MediaQ(vector<double> vettore){
  double somma=0;
  for (int i=0; i<vettore.size(); i++){
    somma += pow(vettore.at(i), 2.);
  }
  return somma/vettore.size();
}

//Varianza dei dati contenuti in un vettore:
double Varianza(vector<double> vettore){
  double somma=0;
  for(unsigned int i=0; i<vettore.size(); i++) {
    somma += pow(Media(vettore) - vettore.at(i), 2.);
  }
  return somma/(vettore.size()-1);
}

//Deviazione standard dei dati:
double Errore(vector<double> media, vector<double> mediaq, int n, int i){
  if (n==0){
    return 0;
  }
  else{
    return sqrt((mediaq.at(n) - pow(media.at(n), 2.))/i);
  }
}

//Funzione per applicare il metodo a blocchi a un vettore e salvarne i risultati:
void Calcolo(int numeroblocchi, vector<double> dati, const char*fileoutput){

	ofstream risultati;
	risultati.open(fileoutput);

	unsigned int numerodati = dati.size();
	int datiperblocco= int(numerodati/numeroblocchi);
  int sottointervalli = datiperblocco/numeroblocchi;
	vector<double> medie(sottointervalli,0);   //Vettore che mi contiene le medie dei singoli blocchi
	vector<double> medieq(sottointervalli,0);   //Vettore che contiene i quadrati delle medie dei singoli blocchi
	vector<double> medieprog(numeroblocchi,0);  //Vettore che mi contiene la media cumulativa dei blocchi
	vector<double> medieqprog(numeroblocchi,0);   //Vettore che contiene i quadrati delle medie cumulative dei blocchi
	vector<double> erroreprog(numeroblocchi,0);      //Vettore che mi contiene gli errori cumulativi

	for(int i=0; i<numeroblocchi; i++){

		for(int j=0; j<sottointervalli; j++){
			double t = sqrt( Somma(dati, i*datiperblocco + j*sottointervalli, i*datiperblocco + (j+1)*sottointervalli)/sottointervalli);

		  medie.at(j) = t;
		  medieq.at(j) = pow(t,2);
	   }

  double mediaprog=0.;
  double mediaqprog=0.;
  double errprog=0.;

  mediaprog = Media(medie);
  mediaqprog = Media(medieq);
  errprog = sqrt((mediaqprog - pow(mediaprog, 2))/(sottointervalli -1));

  if (risultati.is_open()){
		risultati << i << " "<< setprecision(12) <<mediaprog << " " <<setprecision(12) << errprog<< endl;
	 	}
		else cerr << "errore nell'aprire il file dati" << endl;
	}
	risultati.close();
}

//Funzione calcolo chi quadro:
void ChiQuadro(int numintervalli, int numerodati, vector<double> dati, const char* fileoutput){

	ofstream risultati;
  risultati.open(fileoutput);

	//Numero di dati sui cui calcolo il chi quadro singolo
	int datiperintervallo = int(numerodati/numintervalli);

	vector<double> chiquadri(numintervalli,0); //Vettore contenente i chi quadri i-esimi
	vector<double> intervalli(numintervalli+1,0);		//Vettore degli intervalli

	//Vettore che contiene l'estremo destro di ogni intervallo in cui è diviso [0,1)
	for(int i=0; i<numintervalli+1; i++){
		intervalli.at(i) = (1.*i)/numintervalli;
	}

  for(int i=0; i<numintervalli+1; i++){
  //Vettore che conterrà gli osservati per ogni intervallo
  vector<double> osservati(numintervalli,0);

    for(int j=i*datiperintervallo; j<(i+1)*datiperintervallo; j++) {
      for(int s=0; s<numintervalli+1; s++){
        if( (dati.at(j) >= intervalli.at(s) ) && (dati.at(j) < intervalli.at(s+1) ) ){
          osservati.at(s) ++;
        }
      }
    }
		//Calclo del chi quadro per l'i-esimo blocco
    double somma = 0;
    double atteso = ((1.*datiperintervallo)/(1.*numintervalli));		//Conteggi attesi nel signolo blocco
    for(int k=0; k<numintervalli; k++){
      somma +=  pow( osservati.at(k) - atteso ,2) / atteso;
    }
    chiquadri.at(i) = somma;
  }

	for(int i=1; i < chiquadri.size(); i++){
  	if(risultati.is_open()) risultati << i << " " << setprecision(12) << chiquadri.at(i) << endl;
  	else cerr << "Errore nell'aprire il file dei dati" << endl;
  }
  risultati.close();
}
