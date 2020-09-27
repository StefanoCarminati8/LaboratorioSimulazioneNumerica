#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <time.h>
#include "random.h"
#include "metodoblocchi.h"
#include "viaggiatore.h"

using namespace std;

int main (){

	clock_t start,end;
	double time;
	start = clock();

	int numberOfCities = 32;
    int numberOfParents = 1;
    int numberOfMetropolisSteps = 1000;
	int numberOfMutations = 2;

	int numberOfEquilibrationSteps = 1000;
	double InitialBeta = 0.;
	double FinalBeta = 80.;
	double increase = 0.0005;

	TravelingSalesman Walker1(numberOfCities,numberOfParents,numberOfMetropolisSteps);
	double ray = 1.;
	Walker1.Circle(ray);
	double beta = 0.;

	cout << "Città su una circonferenza di raggio 1" << endl;
	Walker1.Equilibration(numberOfEquilibrationSteps, beta, numberOfMutations);
	for(double beta = InitialBeta; beta<FinalBeta; beta += increase){
		Walker1.RunMetropolis(beta, numberOfMutations);
	}
	Walker1.BestPath("path_circle.out");
 	Walker1.BestDistanceAnnealing("distance_circle.out");

	TravelingSalesman Walker2(numberOfCities,numberOfParents,numberOfMetropolisSteps);
	double side = 1.;
	Walker2.Square(side);
	beta = 0.;

	cout << "Città dentro quadrato di lato 1" << endl;
	Walker2.Equilibration(numberOfEquilibrationSteps, beta, numberOfMutations);
	for(double beta = InitialBeta; beta<FinalBeta; beta += increase){
		Walker2.RunMetropolis(beta, numberOfMutations);
	}
	Walker2.BestPath("path_square.out");
  	Walker2.BestDistanceAnnealing("distance_square.out");

	end = clock();
	time = ((double)(end-start))/CLOCKS_PER_SEC;

	cout << endl << "Tempo di esecuzione: " << time << endl;
 	return 0;

}
