#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "metodoblocchi.h"

using namespace std;

double t=0.;
double s0 = 100;
double s=0.;
double T=1.;
double K=100;
double r=0.1;
double sigma=0.25;
double S=0.;
double n=100;
double M=10000;

vector<double> f(M,0);
vector<double> g(M,0);
vector<double> h(M,0);
vector<double> prova(M,0);

double GBM(double t, double T, double f) {
    double S=s0*exp((r-0.5*pow(sigma,2))*(T-t) + f*sqrt(T-t));
    return S;
}

double Call(double t, double T, double N, double f) {
  double step=(T-t)/N;
  for(int i=0; i<N; i++){

    s=GBM(t+i*step,t+(i+1)*step, f);
  //  cout << "valore s dopo: " << s << endl;   //controllo valori
    s0=s;
  //  cout << "valore s0 dopo: " << s0 << endl;
    double call = s*(0.5*(1+erf(((1/(sigma*sqrt(t+(i+1)*step-(t+i*step))))*(log(s/K) + (r+0.5*pow(sigma,2)*(t+(i+1)*step -(t+i*step)))))/sqrt(2))))-K*exp(-r*(t+(i+1)*step-(t+i*step)))*((0.5*(1+erf(((1/(sigma*sqrt(t+(i+1)*step -(t+i*step))))*(log(s/K) + (r+0.5*pow(sigma,2)*(t+(i+1)*step-(t+i*step)))))/sqrt(2))))-sigma*sqrt(t+(i+1)*step-(t+i*step)));
    return call;

  }
}

double Put(double t, double T, double N, double f){
  double step=(T-t)/N;
  for(int i=0; i<N; i++){
    s=GBM(t+i*step,t+(i+1)*step, f);
    s0=s;
    double put = s*((0.5*(1+erf(((1/(sigma*sqrt(t+(i+1)*step-(t+i*step))))*(log(s/K) + (r+0.5*pow(sigma,2)*(t+(i+1)*step -(t+i*step)))))/sqrt(2))))-1)-K*exp(-r*(t+(i+1)*step-(t+i*step)))*(((0.5*(1+erf(((1/(sigma*sqrt(t+(i+1)*step -(t+i*step))))*(log(s/K) + (r+0.5*pow(sigma,2)*(t+(i+1)*step-(t+i*step)))))/sqrt(2))))-sigma*sqrt(t+(i+1)*step-(t+i*step)))-1);
    return put;
  }
}

int main(){

  vector<double> calldirette(M,0);
  vector<double> callastep(M,0);
  vector<double> putdirette(M,0);
  vector<double> putastep(M,0);

  Random rnd;
  SetRandomGenerator(&rnd);
  for(int i =0; i<10000; i++){

    h.at(i) = rnd.Rannyu(0,1);
    calldirette.at(i) = Call(t,T,1,h.at(i));
    callastep.at(i) = Call(t,T,n,h.at(i));
    putdirette.at(i) = Put(t,T,1,h.at(i));
    putastep.at(i) = Put(t,T,n,h.at(i));

    s0=100;
    s=0;

  }

  rnd.SaveSeed();

  Calcolo(n, calldirette, "riscalldirette.out");
  Calcolo(n, callastep, "riscallastep.out");
  Calcolo(n, putdirette, "risputdirette.out");
  Calcolo(n, putastep, "risputastep.out");

  return 0;
}
