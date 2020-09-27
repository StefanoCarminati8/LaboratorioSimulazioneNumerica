#ifndef _viaggiatore_h
#define _viaggiatore_h

#include "random.h"
#include <vector>


using namespace std;

class TravelingSalesman{

private:

  Random _rnd;
  int _NumberOfCities;
  int _NumberOfParents;
  int _NumberOfGenerations;
  int _BestIndex;                             //Indice del parent con la loss function minore
  int _StepsAccepted;                         //Passi del metropolis accettati

  vector<vector<int>> _StringsOfCities;       //Matrice delle posizioni della città dell'intera popolazione
  vector<vector<double>> _Positions;          //Matrice corrispondente di tutte le posizioni delle città
  vector<double> _FitnessDistances;           //Vettore che contiene i valori delle loss functions per tutti i parents
  vector <int> _Parents;                      //Vettori degli indici dei parents scelti dal metodo di selezione
  vector<vector<int>> _NewGeneration;         //Matrice della prole
  vector<double> _RememberBestDistance;       //Vettore che tiene in imemoria la minor distanza per ogni generazione
  vector<vector<double>> _BestDistance50;     //Matrice della media e std delle migliori distanze ottenute da 50% popolazione

public:

  TravelingSalesman(int NumberOfCity, int NumberOfParents, int NumberOfGenerations);
  ~TravelingSalesman();

  void Restart(void);

  //Generazione delle città
  void Square(double ray);
  void Circle(double side);

  void Fitness(void);
  //Selezioni
  void ParentsSelection(void);
  int Selection(void);
  void ParentsOrderSelection(void);
  int OrderSelection(void);
  //CrossOver
  void OneOffspringCrossover(int parent1, int parent2);
  void TwoOffspringCrossover(int parent1, int parent2, double ProbabilityOfCrossover);
  void NoCrossover(void);
  //Mutazioni
  void MutationPairPermutation(double ProbabilityOfMutation);
  void MutationTotCitiesPermutation(double ProbabilityOfMutation);
  void MutationShiftAll(double ProbabilityOfMutation);
  void MutationShiftTotCities(double ProbabilityOfMutation);
  void SelectMutation(double ProbabilityOfMutation);
  //CheckFunction
  void CheckBounds(void);

  void RunGA(double ProbabilityOfMutation, double ProbabilityOfCrossover);

  //Risultati
  void Print(void);                             //Stampa a schermo le citta per tutti i genitori
  void PrintGeneration(void);
  void BestPath(const char * outputFile);
  void BestDistance(const char * outputFile);
  void UpdateDistance50(int generation);
  void BestDistance50(const char * outputFile);
  void BestIndexFitness(void);

  //Simulated Annealing
  vector<double> FitnessAnnealing(vector<vector<int>> matrix);

  void MutationPairPermutationAnnealing(int NumberOfMutations);
  void MutationTotCitiesPermutationAnnealing(int NumberOfMutations);
  void MutationShiftAllAnnealing(int NumberOfMutations);
  void MutationShiftTotCitiesAnnealing(int NumberOfMutations);
  void SelectMutationAnnealing(int NumberOfMutations);

  void PrintAcceptanceRate();
  void Equilibration(int NumberOfEquilibrationSteps, double Beta, int NumberOfMutations);
  void Metropolis(double Beta, int NumberOfMutations);
  void RunMetropolis(double Beta, int NumberOfMutations);
  void BestDistanceAnnealing(const char * outputFile);

};

#endif
