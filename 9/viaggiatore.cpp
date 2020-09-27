#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>              //random_shuffle
#include "random.h"
#include "metodoblocchi.h"
#include "viaggiatore.h"


using namespace std;

//Funzione che restituisce il segno del numero
double Signum (double x){
  if (x>=0){return 1.;} else{return -1.;}
}

//Costruttore
TravelingSalesman::TravelingSalesman(int NumberOfCities, int NumberOfParents, int NumberOfGenerations): _FitnessDistances(NumberOfParents){

  SetRandomGenerator(&_rnd);
  _NumberOfGenerations= NumberOfGenerations;
  _NumberOfCities = NumberOfCities;
  _NumberOfParents = NumberOfParents;

  _BestIndex = 0;

  //Configurazione standard [0,1,2,......,29]
  vector<int> cities(_NumberOfCities);
  for (int i=0; i<_NumberOfCities; i++) cities.at(i) = i;

  //Ottengo la popolazione iniziale attraverso uno shuffle della configurazione standard con il constraint che la prima città e l'ultima devono essere sempre la stessa
  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());   //random_shuflle[beginIndex, lastIndex)
    _StringsOfCities.push_back(cities);             //Ottengo la matrice delle posizioni di tutte le città
  }
  _NewGeneration = _StringsOfCities;
  _StepsAccepted = 0;
}

//Distruttore
TravelingSalesman::~TravelingSalesman(){
  _rnd.SaveSeed();
}

//Check function per verificare i constraints
void TravelingSalesman::CheckBounds(void){

  for(int i=0; i<_NumberOfParents; i++){
    vector<int> auxiliary;
    for(int j=0; j<_NumberOfCities; j++){
        auxiliary.push_back(_StringsOfCities[i][j]);
    }
    vector<int>::iterator ip;
    ip = unique(auxiliary.begin(), auxiliary.end());
    auxiliary.resize(std::distance(auxiliary.begin(), ip));
    if(auxiliary.size() != _NumberOfCities) {
      cout << endl << "PROBLEM: the parents no. " << i << " does not have all the different cities" << endl;
      cout << endl;
      cout << "Parent with problems in detail: ";
      for (int j=0; j<_NumberOfCities; j++){
        cout << _StringsOfCities[i][j] << " ";
      }
      cout << endl;
    }
  }
}

//Metodo che valuta la loss function per tutta la generazione
void TravelingSalesman::Fitness(void){
  for(int j=0; j<_NumberOfParents; j++){
    double lossFunction = 0;
    for (int i=0; i<_NumberOfCities-1;i++){
      int index = _StringsOfCities[j][i];
      int indexPlus = _StringsOfCities[j][i+1];
      lossFunction += pow(_Positions[index][1]-_Positions[indexPlus][1],2) + pow(_Positions[index][2]-_Positions[indexPlus][2],2);
    }
    //Aggiungo la distanze per tornare alla hometown ----> è come se la 31-esima città fosse la hometown (soddisfo il constraints non testato)
    lossFunction += pow(_Positions[_StringsOfCities[j][0]][1]-_Positions[_StringsOfCities[j][_NumberOfCities-1]][1],2) + pow(_Positions[_StringsOfCities[j][0]][2]-_Positions[_StringsOfCities[j][_NumberOfCities-1]][2],2);
    _FitnessDistances[j] = lossFunction;
  }
}

//Metodo che seleziona tra tutti gli individui i possibili genitori
void TravelingSalesman::ParentsSelection(void){

  vector <int> parents;
  vector <double> orderedDistances = _FitnessDistances;
  sort(orderedDistances.begin(),orderedDistances.end());

  for(int i=0;i<int(orderedDistances.size()/4); i++){
    double possible = orderedDistances.at(i);
    vector<double>::iterator indexParent = find(_FitnessDistances.begin(), _FitnessDistances.end(), possible);
    parents.push_back(distance(_FitnessDistances.begin(), indexParent));
  }
  for(int i=int(orderedDistances.size()/5); i<int(orderedDistances.size()/2); i++){
    double random = int(_rnd.Rannyu(orderedDistances.size()/5,orderedDistances.size()/2));
    double possible = orderedDistances.at(random);
    vector<double>::iterator indexParent = find(_FitnessDistances.begin(), _FitnessDistances.end(), possible);
    parents.push_back(distance(_FitnessDistances.begin(), indexParent));
  }
  _Parents=parents;
}

//Metodo che seleziona tra i genitori possibili uno casuale
int TravelingSalesman::Selection(void){
  double p = int(_rnd.Rannyu(0,_Parents.size()));
  return _Parents.at(p);
}

//Metodo che seleziona tra tutti gli individui i possibili genitori
void TravelingSalesman::ParentsOrderSelection(void){

  vector <int> parents;
  vector <double> orderedDistances = _FitnessDistances;
  sort(orderedDistances.begin(),orderedDistances.end());

  for(int i=0;i<int(orderedDistances.size());i++){
    double possible = orderedDistances.at(i);
    vector<double>::iterator indexParent = find(_FitnessDistances.begin(), _FitnessDistances.end(), possible);
    parents.push_back(distance(_FitnessDistances.begin(), indexParent));
  }
  _Parents = parents;
}
//Metodo che fa una selezione in base alla migliore distanza
int TravelingSalesman::OrderSelection(void){

  double p = _rnd.Rannyu();
  return _Parents.at(_NumberOfCities * int(pow(p,5)));

}

void TravelingSalesman::OneOffspringCrossover(int parent1, int parent2){

  unsigned int indexStart = int(_rnd.Rannyu(0,_NumberOfCities));
  unsigned int indexEnd = int(_rnd.Rannyu(indexStart,_NumberOfCities));

  vector<int> offspring(_NumberOfCities, 0);
  //Prima parte della prole
  for(unsigned int i=indexStart; i<indexEnd; i++) offspring.at(i) = _StringsOfCities[parent1][i];
  int actualPosition = 0;
  for(unsigned int i=0; i<indexStart; i++){
    for(unsigned int j=actualPosition; j<_NumberOfCities; j++){
      int counter = 0;
      for(unsigned int k=indexStart; k<indexEnd; k++){
        if(_StringsOfCities[parent2][j] == offspring.at(k)) counter++;
      }
      if(counter == 0){
        offspring.at(i) = _StringsOfCities[parent2][j];
        actualPosition = j+1;
        break;
      }
    }
  }

  //Seconda parte del figlio
  for(int i=indexEnd; i<_NumberOfCities; i++){
    for(unsigned int j=actualPosition; j<_NumberOfCities; j++){
      int counter = 0;
      for(unsigned int k=indexStart; k<indexEnd; k++){
        if(_StringsOfCities[parent2][j] == offspring.at(k)) counter++;
      }
      if(counter == 0){
        offspring.at(i) = _StringsOfCities[parent2][j];
        actualPosition = j+1;
        break;
      }
    }
  }
  _NewGeneration.push_back(offspring);

}

//Metodo diverso dal precedente crossover che avviene con una certa probabilità impostata esternamente e ottiene due figli
void TravelingSalesman::TwoOffspringCrossover(int parent1, int parent2, double probabilityOfCrossover){

  unsigned int indexStart = int(_rnd.Rannyu(1,_NumberOfCities));
  vector<int> offspring1(_NumberOfCities, 0);
  vector<int> offspring2(_NumberOfCities, 0);

  double random = _rnd.Rannyu();
  if(random<probabilityOfCrossover){
    for(unsigned int i=0; i<indexStart; i++) offspring1.at(i) = _StringsOfCities[parent1][i];
    for(unsigned int i=0; i<indexStart; i++) offspring2.at(i) = _StringsOfCities[parent2][i];
    //Primo figlio
    int actualPosition = 0;
    for(unsigned int i=indexStart; i<_NumberOfCities; i++){
      for(unsigned int j=actualPosition; j<_NumberOfCities; j++){
        int counter = 0;
        for(unsigned int k=0; k<indexStart; k++){
          if(_StringsOfCities[parent2][j] == offspring1.at(k)) counter++;
        }
        if(counter == 0){
          offspring1.at(i) = _StringsOfCities[parent2][j];
          actualPosition = j+1;
          break;
        }
      }
    }
    _NewGeneration.push_back(offspring1);

    //Secondo figlio
    actualPosition = 0;
    for(unsigned int i=indexStart; i<_NumberOfCities; i++){
      for(unsigned int j=actualPosition; j<_NumberOfCities; j++){
        int counter = 0;
        for(unsigned int k=0; k<indexStart; k++){
          if(_StringsOfCities[parent1][j] == offspring2.at(k)) counter++;
        }
        if(counter == 0){
          offspring2.at(i) = _StringsOfCities[parent1][j];
          actualPosition = j+1;
          break;
        }
      }
    }
    _NewGeneration.push_back(offspring2);
  }
  else _NewGeneration = _StringsOfCities;
}

//Metodo per testare il GA senza crossover
void TravelingSalesman::NoCrossover(void){
  _NewGeneration = _StringsOfCities;
}

//Metodo che permuta una coppia di città all'interno di tutti i parent
void TravelingSalesman::MutationPairPermutation(double ProbabilityOfMutation){
  double random = _rnd.Rannyu();
  if(random<ProbabilityOfMutation){
    for(unsigned int i=0; i <_NewGeneration.size();i++){
      unsigned int indexStart = int(_rnd.Rannyu(0,_NumberOfCities));
      unsigned int indexEnd = int(_rnd.Rannyu(0,_NumberOfCities));
      int auxiliary = _NewGeneration[i][indexStart];
      _NewGeneration[i][indexStart] = _NewGeneration[i][indexEnd];
      _NewGeneration[i][indexEnd] = auxiliary;
    }
  }
}

//Metodo che permuta un gruppo di città all'interno di tutti i parent
void TravelingSalesman::MutationTotCitiesPermutation(double ProbabilityOfMutation){

  double random = _rnd.Rannyu();
  if(random<ProbabilityOfMutation){
    for(unsigned int i=0; i <_NewGeneration.size();i++){
      vector<int> auxiliary;
      double numberOfCities = int(_rnd.Rannyu(1, int(_NumberOfCities/2)));
      for(unsigned int j=0; j<numberOfCities; j++) auxiliary.push_back(_NewGeneration[i][j]);
      for(unsigned int j=0; j<numberOfCities; j++) _NewGeneration[i][j] = _NewGeneration[i][_NumberOfCities - (numberOfCities - j)];
      for(unsigned int j=0; j<numberOfCities; j++) _NewGeneration[i][_NumberOfCities - (numberOfCities - j)] = auxiliary[j];
    }
  }
}

//Metodo che fa uno shift di tutte le città di un certo valore scelto casualmente
void TravelingSalesman::MutationShiftAll(double ProbabilityOfMutation){

  vector<vector <int>> shiftGeneration;
  shiftGeneration = _NewGeneration;
  int NumberOfShift = int(_rnd.Rannyu(0,_NumberOfCities));
  double random = _rnd.Rannyu();
  if(random<ProbabilityOfMutation){
    for(unsigned int i=0; i<_NumberOfParents;i++){
      for(unsigned int j=0; j<_NumberOfCities;j++){
        shiftGeneration[i][(j+NumberOfShift)%_NumberOfCities] = _NewGeneration[i][j];
      }
    }
    _NewGeneration = shiftGeneration;
  }
}

//Metodo che shifta un numero di città ogni volta diverso  di un certo valore scelto casualmente
void TravelingSalesman::MutationShiftTotCities(double ProbabilityOfMutation){

  //Escludo la hometown e quelle per cui necessitano le PBC
  int NumberOfShift = int(_rnd.Rannyu(0,_NumberOfCities));
  unsigned int indexStart = int(_rnd.Rannyu(1,_NumberOfCities - (NumberOfShift+1)));
  vector<vector <int>> shiftGeneration;
  shiftGeneration = _NewGeneration;
  double random = _rnd.Rannyu();
  if(random<ProbabilityOfMutation){
    for(unsigned int i=0; i<_NumberOfParents;i++){
      for(unsigned int j=indexStart; j<_NumberOfCities;j++){
        if(j<(_NumberOfCities - NumberOfShift))  shiftGeneration[i][(j+NumberOfShift)%(_NumberOfCities)] = _NewGeneration[i][j];
        else shiftGeneration[i][(j+ NumberOfShift + indexStart)%(_NumberOfCities)] = _NewGeneration[i][j];
      }
    }
    _NewGeneration = shiftGeneration;
  }

}

//Metodo che scegli casualmente quale mutazione applicare
void TravelingSalesman::SelectMutation(double ProbabilityOfMutation){

  int mutation = _rnd.Rannyu(0,4);
  if (mutation == 0) MutationPairPermutation(ProbabilityOfMutation);
  if (mutation == 1) MutationTotCitiesPermutation(ProbabilityOfMutation);
  if (mutation == 2) MutationShiftAll(ProbabilityOfMutation);
  if (mutation == 3) MutationShiftTotCities(ProbabilityOfMutation);

}

//Run di tutte le generazioni
void TravelingSalesman::RunGA(double ProbabilityOfMutation, double ProbabilityOfCrossover){

  for(int i=0; i <_NumberOfGenerations;i++){
    if(i%50 == 0) cout << "GENERAZIONE NUMERO: " << i << endl;
    BestIndexFitness();
    UpdateDistance50(i);
    _RememberBestDistance.push_back(_FitnessDistances.at(_BestIndex));

    //Test dell GA senza crossover ma solo con mutazioni
    //NoCrossover();

    //Crossover restituendo un solo figlio dati due genitori
    for(int j=0; j<_NumberOfParents; j++){
    ParentsSelection();
    int parent1 = Selection();
    int parent2 = Selection();
    OneOffspringCrossover(parent1, parent2);
    }
/*
    //Crossover restituendo due figli dati due genitori
    for(int j=0; j<int(_NumberOfParents/2); j++){
      ParentsSelection();
      int parent1 = Selection();
      int parent2 = Selection();
      TwoOffspringCrossover(parent1, parent2, ProbabilityOfCrossover);
    }
*/
    //MUTAZIONE

    //Scelgo random una mutazione da fare ad ogni generazione

    //SelectMutation(ProbabilityOfMutation);

    //Mutazioni singole
    MutationPairPermutation(ProbabilityOfMutation);
    //MutationTotCitiesPermutation(ProbabilityOfMutation);
    //MutationShiftAll(ProbabilityOfMutation);
    //MutationShiftTotCities(ProbabilityOfMutation);

    CheckBounds();
    _StringsOfCities = _NewGeneration;
    _NewGeneration.erase(_NewGeneration.begin(),_NewGeneration.end());
  }
  BestIndexFitness();
  return;
}

//Metodo che mi identifica l'indice del parent con la loss function minore
void TravelingSalesman::BestIndexFitness(void){
  Fitness();
  _BestIndex = min_element(_FitnessDistances.begin(), _FitnessDistances.end()) - _FitnessDistances.begin();
}

//Metodo che stampa in un file il cammino con la loss function minore
void TravelingSalesman::BestPath(const char * outputFIle){

  ofstream results;
  results.open(outputFIle);

  BestIndexFitness();
  //Stampo  x, y delle città in ordine di passaggio dalla prima all'ultima
  for (int i=0; i <_NumberOfCities;i++){
    results << _Positions[_StringsOfCities[_BestIndex][i]][1] << " " << _Positions[_StringsOfCities[_BestIndex][i]][2]  << endl;
  }
  //Stampo come posizione finale la hometown
  results << _Positions[_StringsOfCities[_BestIndex][0]][1] << " " << _Positions[_StringsOfCities[_BestIndex][0]][2] <<  endl;
  results.close();
}

//Metodo che stampa in un file la distanza minore per ogni generazione
void TravelingSalesman::BestDistance(const char * outputFIle){
  ofstream results;
  results.open(outputFIle);

  for(int i=0; i<_NumberOfGenerations;i++){
    results <<  _RememberBestDistance[i] << endl;
  }
  results.close();
}

//Metodo che calcola media e std delle migliori distanze del 50% popolazione
void TravelingSalesman::UpdateDistance50(int generation) {
  vector <double> orderedDistances = _FitnessDistances;
  sort(orderedDistances.begin(),orderedDistances.end());
  vector<double> distances50(int(_NumberOfParents/2));
  for(int i = 0; i< int(_NumberOfParents/2); i++){
    distances50[i] = orderedDistances[i];
  }
  vector<double> auxiliary(3);
  auxiliary[0] = generation;
  auxiliary[1] = Mean(distances50);
  auxiliary[2] = sqrt(ErrorVector(distances50));
  _BestDistance50.push_back(auxiliary);
  return;
}

//Metodo che scrive su file le migliori distanze del 50% popolazione di ogni generazione
void TravelingSalesman::BestDistance50(const char * outputFile) {
  ofstream results;
  results.open(outputFile);

  for(unsigned int i=0; i<_BestDistance50.size();i++){
    results << _BestDistance50[i][0] << " " << _BestDistance50[i][1] << " " << _BestDistance50[i][2] << endl;
  }
  results.close();
  return;

}

//Metodo  che mi riinizializza tutto e mi genera una nuova generazione iniziale
void TravelingSalesman::Restart(void){
  _StringsOfCities.erase(_StringsOfCities.begin(),_StringsOfCities.end());
  _NewGeneration.erase(_NewGeneration.begin(),_NewGeneration.end());
  _Positions.erase(_Positions.begin(),_Positions.end());
  _BestIndex = 0;
  _RememberBestDistance.erase(_RememberBestDistance.begin(),_RememberBestDistance.end());
  srand(10);
  vector<int> cities(_NumberOfCities);
  for (int i=0; i<_NumberOfCities; i++){
    cities.at(i) = i;
  }

  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());
    _StringsOfCities.push_back(cities);
  }
}

//Generazione delle città
//Città all'interno di un quadrato
void TravelingSalesman::Square(double side){
  vector<double> positionsOfCities(3);
  for (int i=0; i<_NumberOfCities;i++){

    double x = _rnd.Rannyu(0,side);
    double y = _rnd.Rannyu(0,side);

    positionsOfCities[0] = i; positionsOfCities[1] = x; positionsOfCities[2] = y;     //Numero città ---- Coordinata x ---- Coordinata y
    _Positions.push_back(positionsOfCities);
  }

}

//Città su una circonferenza
void TravelingSalesman::Circle(double ray){
  vector<double> positionsOfCities(3);
  for (int i=0; i<_NumberOfCities;i++){

    double x = _rnd.Rannyu(-ray,ray);
    double y = sqrt(ray*ray -x*x) * Signum(_rnd.Rannyu(-1,1));

    positionsOfCities[0] = i; positionsOfCities[1] = x; positionsOfCities[2] = y;     //Numero città ---- Coordinata x ---- Coordinata y
    _Positions.push_back(positionsOfCities);
  }
}

//Metodo per stampare le posizioni delle città
void TravelingSalesman::Print(void){

  for(int j=0; j<_NumberOfParents; j++){
		for(int i=0; i<_NumberOfCities;i++){
			cout << _StringsOfCities[j][i] << " ";
		}
    cout << endl;
	}
}

//Metodo per stampare la nuova generazione
void TravelingSalesman::PrintGeneration(void){

  for(int j=0; j<_NumberOfParents; j++){
		for(int i=0; i<_NumberOfCities;i++){
			cout << _NewGeneration[j][i] << " ";
		}
    cout << endl;
	}

}

vector<double> TravelingSalesman::FitnessAnnealing(vector<vector<int>> matrix){

  vector<double> distances(_NumberOfParents);
  for(int j=0; j<_NumberOfParents; j++){
    double lossFunction = 0;
    for (int i=0; i<_NumberOfCities-1;i++){
      int index = matrix[j][i];
      int index_plus = matrix[j][i+1];
      lossFunction += pow(_Positions[index][1]-_Positions[index_plus][1],2) + pow(_Positions[index][2]-_Positions[index_plus][2],2);
    }
    lossFunction += pow(_Positions[matrix[j][0]][1]-_Positions[matrix[j][_NumberOfCities-1]][1],2) + pow(_Positions[matrix[j][0]][2]-_Positions[matrix[j][_NumberOfCities-1]][2],2);
    distances[j] = lossFunction;
  }
  return distances;
}


//Metodo che permuta una coppia di città all'interno di tutti i parent
void TravelingSalesman::MutationPairPermutationAnnealing(int NumberOfMutations){

  for(int j=0; j<NumberOfMutations; j++){
    for(unsigned int i=0; i <_NewGeneration.size();i++){
      unsigned int indexStart = int(_rnd.Rannyu(0,_NumberOfCities));
      unsigned int indexEnd = int(_rnd.Rannyu(0,_NumberOfCities));
      int auxiliary = _NewGeneration[i][indexStart];
      _NewGeneration[i][indexStart] = _NewGeneration[i][indexEnd];
      _NewGeneration[i][indexEnd] = auxiliary;
    }
  }
}

//Metodo che permuta un gruppo di città all'interno di tutti i parent
void TravelingSalesman::MutationTotCitiesPermutationAnnealing(int NumberOfMutations){

  for(int j=0; j<NumberOfMutations; j++){
    for(unsigned int i=0; i <_NewGeneration.size();i++){
      vector<int> auxiliary;
      double numberOfCities = int(_rnd.Rannyu(1, int(_NumberOfCities/2)));
      for(unsigned int j=0; j<numberOfCities; j++) auxiliary.push_back(_NewGeneration[i][j]);
      for(unsigned int j=0; j<numberOfCities; j++) _NewGeneration[i][j] = _NewGeneration[i][_NumberOfCities - (numberOfCities - j)];
      for(unsigned int j=0; j<numberOfCities; j++) _NewGeneration[i][_NumberOfCities - (numberOfCities - j)] = auxiliary[j];
    }
  }
}

//Metodo che fa uno shift di tutte le città di un certo valore scelto casualmente
void TravelingSalesman::MutationShiftAllAnnealing(int NumberOfMutations){

  vector<vector <int>> shiftGeneration;
  shiftGeneration = _NewGeneration;
  int NumberOfShift = int(_rnd.Rannyu(0,_NumberOfCities));
  for(int j=0; j<NumberOfMutations; j++){
    for(unsigned int i=0; i<_NumberOfParents;i++){
      for(unsigned int j=0; j<_NumberOfCities;j++){
        shiftGeneration[i][(j+NumberOfShift)%_NumberOfCities] = _NewGeneration[i][j];
      }
    }
    _NewGeneration = shiftGeneration;
  }
}

//Metodo che shifta un numero di città ogni volta diverso  di un certo valore scelto casualmente
void TravelingSalesman::MutationShiftTotCitiesAnnealing(int NumberOfMutations){

  //Escludo la hometown e quelle per cui necessitano le PBC
  int NumberOfShift = int(_rnd.Rannyu(0,_NumberOfCities));
  unsigned int indexStart = int(_rnd.Rannyu(1,_NumberOfCities - (NumberOfShift+1)));
  vector<vector <int>> shiftGeneration;
  shiftGeneration = _NewGeneration;
  for(int j=0; j<NumberOfMutations; j++){
    for(unsigned int i=0; i<_NumberOfParents;i++){
      for(unsigned int j=indexStart; j<_NumberOfCities;j++){
        if(j<(_NumberOfCities - NumberOfShift))  shiftGeneration[i][(j+NumberOfShift)%(_NumberOfCities)] = _NewGeneration[i][j];
        else shiftGeneration[i][(j+ NumberOfShift + indexStart)%(_NumberOfCities)] = _NewGeneration[i][j];
      }
    }
    _NewGeneration = shiftGeneration;
  }

}

//Metodo che scegli casualmente quale mutazione applicare
void TravelingSalesman::SelectMutationAnnealing(int NumberOfMutations){

  int mutation = _rnd.Rannyu(0,4);
  if (mutation == 0) MutationPairPermutationAnnealing(NumberOfMutations);
  if (mutation == 1) MutationTotCitiesPermutationAnnealing(NumberOfMutations);
  if (mutation == 2) MutationShiftAllAnnealing(NumberOfMutations);
  if (mutation == 3) MutationShiftTotCitiesAnnealing(NumberOfMutations);

}

void TravelingSalesman::Metropolis(double Beta, int NumberOfMutations){

  BestIndexFitness();

  _RememberBestDistance.push_back(_FitnessDistances.at(_BestIndex));

  vector<double> oldFitness = FitnessAnnealing(_StringsOfCities);

  //Mutazioni singole
  //MutationPairPermutationAnnealing(NumberOfMutations);
  //MutationTotCitiesPermutation(NumberOfMutations);
  //MutationShiftAllAnnealing(NumberOfMutations);
  //MutationShiftTotCitiesAnnealing(NumberOfMutations);

  //Scelta random di una mutazione
  SelectMutationAnnealing(NumberOfMutations);

  vector<double> newFitness = FitnessAnnealing(_NewGeneration);

  for(unsigned int i=0; i<_FitnessDistances.size(); i++){
    double boltzmannWeight = exp(-Beta * (newFitness[i] - oldFitness[i]));
    if(boltzmannWeight >= 1){
	    _StringsOfCities[i] = _NewGeneration[i];
	    _StepsAccepted++;
    }
    else {
	    double random = _rnd.Rannyu();
	    if(random <= boltzmannWeight){
		    _StringsOfCities[i] = _NewGeneration[i];
	  	  _StepsAccepted++;
	    }
    }
  }
  _NewGeneration = _StringsOfCities;
  return;
}

void TravelingSalesman::RunMetropolis(double Beta, int NumberOfMutations) {

  BestIndexFitness();
	for(int i=0; i<_NumberOfGenerations; i++) {
		Metropolis(Beta, NumberOfMutations);
	}
	BestIndexFitness();
}

void TravelingSalesman::Equilibration(int NumberOfEquilibrationSteps, double Beta, int NumberOfMutations) {

	for(int i=0; i<NumberOfEquilibrationSteps; i++) {
		Metropolis(Beta, NumberOfMutations);
	}
	_RememberBestDistance.erase(_RememberBestDistance.begin(),_RememberBestDistance.end());
	_StepsAccepted = 0;
}

void TravelingSalesman::PrintAcceptanceRate(){
	cout << "Acceptance rate: " << ((((double)_StepsAccepted/(double)_NumberOfGenerations))/(double)_NumberOfParents)*100. << "%" << endl;
  _StepsAccepted = 0;
}

void TravelingSalesman::BestDistanceAnnealing(const char * outputFile){
  ofstream results;
  results.open(outputFile);
  vector<double> auxiliary(_NumberOfGenerations);

  for(int j=0; j<_RememberBestDistance.size()/_NumberOfGenerations; j++ ){
    for(int i=0; i< _NumberOfGenerations; i++){
      auxiliary[i] = _RememberBestDistance[j*_NumberOfGenerations+i];
    }
    int minimum = min_element(auxiliary.begin(),auxiliary.end()) - auxiliary.begin();
    results << auxiliary[minimum]<< endl;
  }
    results.close();
}
