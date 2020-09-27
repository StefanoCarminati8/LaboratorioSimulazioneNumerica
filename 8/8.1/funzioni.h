#ifndef _funzioni_h
#define _funzioni_h
#include <cmath>
#include <vector>

using namespace std;

// Classe contentente tutte le funzioni che dovrò intergrare

class Function {
	public:
		virtual double Eval(vector <double> pos) const=0;
		virtual ~Function();
};


// Distribuzioni di probabilità da campionare
class WaveFunction: public Function {

	private:
		double _mi;
		double _sigma;

	public:

		WaveFunction();
		WaveFunction(double mi, double sigma);
		~WaveFunction();

		virtual double Eval(vector <double> pos) const;
		void SetMi(double mi) { _mi = mi; };
		double GetMi() {return _mi;};
		void SetSigma(double sigma) { _sigma = sigma; };
		double GetSigma() {return _sigma;};
	};

class Hamiltonian: public Function{

	private:
		WaveFunction *_wavefunction;

	public:

		Hamiltonian(WaveFunction *wavefunction);
		~Hamiltonian();
		virtual double Eval(vector <double> pos) const;

	};





#endif
