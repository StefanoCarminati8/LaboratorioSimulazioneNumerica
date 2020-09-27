#ifndef _funzioni_h
#define _funzioni_h
#include <cmath>

// Classe contentente tutte le funzioni che dovrò intergrare

class Funzione {
	public:
		virtual double Eval(double x) const=0;
		virtual ~Funzione() = 0;
};


// Funzione da integrare nel primo punto dell' esercizio 02.1
class Coseno: public Funzione {

	public:
		virtual double Eval(double x) const;
		~Coseno();
};

//Funzione da integrare nel secondo punto dell'esercizio 02.1
class gfunz: public Funzione {

	public:
		virtual double Eval(double x) const;
		~gfunz();
};

#endif
