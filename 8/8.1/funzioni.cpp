#include <cmath>
#include <vector>
#include <iostream>
#include "funzioni.h"

using namespace std;

Function::~Function(){}

//Wavefunction
WaveFunction :: WaveFunction(){

	_mi = 0;
	_sigma = 0;

}


WaveFunction :: WaveFunction(double mi, double sigma){

	_mi = mi;
	_sigma = sigma;

}

WaveFunction::~WaveFunction(){}

double WaveFunction :: Eval(vector <double> pos) const {

	double exp_dx = exp(-pow(pos.at(0)-_mi,2) / (2*_sigma*_sigma));
	double exp_sx = exp(-pow(pos.at(0)+_mi,2) / (2*_sigma*_sigma));

	return exp_sx + exp_dx;

}

//Hamiltonian
Hamiltonian :: Hamiltonian(WaveFunction *wavefunction)	{
	_wavefunction = wavefunction;
}

Hamiltonian::~Hamiltonian(){};


double Hamiltonian :: Eval(vector <double> pos) const{		//(H applicato a psi) / psi

	double exp_dx = exp(-pow(pos.at(0)-_wavefunction->GetMi(),2) / (2*_wavefunction->GetSigma()*_wavefunction->GetSigma()));
	double exp_sx = exp(-pow(pos.at(0)+_wavefunction->GetMi(),2) / (2*_wavefunction->GetSigma()*_wavefunction->GetSigma()));

	//Variabili di appoggio
	double fac_1 = 1 / (_wavefunction->GetSigma()*_wavefunction->GetSigma());
	double fac_2 = 1 / pow(_wavefunction->GetSigma(),4);
	double arg_1 = pow(pos.at(0)-_wavefunction->GetMi(),2);
	double arg_2 = pow(pos.at(0)+_wavefunction->GetMi(),2);

	double kineticEnergy = -0.5*(-exp_dx*fac_1 - exp_sx*fac_1 + exp_dx*fac_2*arg_1 + exp_sx*fac_2*arg_2)/_wavefunction->Eval(pos);
	double potential = pow(pos.at(0),4) - 2.5*pow(pos.at(0),2);

	return kineticEnergy + potential;

}
