#ifndef _funzioni_h
#define _funzioni_h
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

class Function {
	public:
		virtual double Eval(vector <double> pos) const=0;
		virtual ~Function();
};

// Distribuzioni di probabilità da campionare
class Uno_s: public Function {
	public:
		virtual double Eval(vector <double> pos) const;
		~Uno_s();
};

class Due_p: public Function {

	public:
		virtual double Eval(vector <double> pos) const;
		~Due_p();
};

#endif
