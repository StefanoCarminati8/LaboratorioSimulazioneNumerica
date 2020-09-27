#ifndef _Integrale_h
#define _Integrale_h
#include "random.h"
#include "funzioni.h"

class Integrale{

  private:

    Random _rnd;
    //Utile per tenere conto di eventuali scambi negli estremi di integrazione
    int _Sign;
    //Estremi di integrazione
    double _a, _b;
    double _sommaparziale;
    double _Integrale;
    Funzione * _integranda;

  public:

    //Costruttore
    Integrale(double a, double b, Funzione * Funzione);
    //Distruttore
  	~Integrale();
    //Calcola l'integrale montecarlo generando numeri uniformenete distribuiti tra _a e _b suddiviso in NumberOfIntervals intervalli con in metodo della media
  	double Media(int NumberOfIntervals);
    // Calcola l'integrale con numeri uniformemente distribuiti tra 0 e 1 con p(x) = 2(1-x)
    double MediaRetta(int NumberOfIntervals);
    // Restituisce il risultato salvato in precedenza
  	double GetResult() const;
    void Start();

};

#endif
