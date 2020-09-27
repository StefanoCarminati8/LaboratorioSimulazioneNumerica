#ifndef _metodoblocchi_h_
#define _metodoblocchi_h_

#include "random.h"
#include <vector>

using namespace std;

//Somma di elementi di un vettore:
double Somma(vector<double> vettore, int inizio, int fine);
//Media dei dati contenuti in un vettore:
double Media(vector<double> vettore);
//Media dei quadrati dei dati contenuti in un vettore:
double MediaQ(vector<double> vettore);
//Varianza dei dati contenuti in un vettore:
double Varianza(vector<double> vettore);
//Deviazione standard dei dati:
double Errore(vector<double> media, vector<double> mediaq, int n, int i);
//Funzione per applicare il metodo a blocchi a un vettore e salvarne i risultati:
void Calcolo(int numerblocchi, vector<double> dati, const char*fileoutput);
//Calcolo Chi Quadro:
void ChiQuadro(int numerodati, int numeroblocchi, vector<double> dati, const char*fileoutput);
#endif
