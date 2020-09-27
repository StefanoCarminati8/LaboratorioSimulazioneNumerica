#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input(); //Inizialization
  int nstepeq = 1000;
  cout << "Equilibrazione per le misure:..." << endl;
  Equilibrazione(nstepeq);
    do{
      for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
      //  accepted=0;
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move(metro);
          Measure();
          Accumulate(); //Update block averages
        }
        if(iblk == nblk){
          cout << "temperatura: " << temp << endl;
        }
        if(metro==1){
      //    cout << "rate accettazione: " << accepted/attempted << endl << endl;
        }
        Averages(iblk);
      }
      //ConfFinal(); //Write final configuration
      temp += 0.05;
    }
    while(temp <= 2.05);

  cout << "-------- misure con campo esterno ------------ " << endl;

  temp = 0.5;
  h=0.02;
  beta=1/temp;
  cout << "equilibrazione del sistema..." << endl;
  Equilibrazione(nstepeq);
    do{
      beta = 1/temp;
      for(int iblk=1; iblk<=nblk; iblk++){
        Reset(iblk);
        for(int istep=1; istep<=nstep; istep++){
          Move(metro);
          Measure();
          Accumulate();
        }
        if(iblk == nblk){
          cout << "temperatura: " << temp << endl;
          if(metro==1){
      //      cout << "rate accettazione: " << accepted/attempted << endl << endl;
          }
        }
        Magnetizzazzione2(iblk);
        if(iblk == nblk){
          cout << "--------------" << endl;
        }
      }
      temp += 0.05;
    }
  while(temp <= 2.05);
  ConfFinal();
  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> caso;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

//initial configuration

if(caso == 1) {
  cout << "parto da configurazione old.0" << endl << endl;
  ifstream ReadConfig;
  ReadConfig.open("old.0");

  for(int i=0; i<nspin; i++){
    ReadConfig >> s[i];
  }
  ReadConfig.close();
}

else if(caso == 0);
cout << "configurazione iniziale casuale: " << endl << endl;
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      Metropolis(o);

    }
    else //Gibbs sampling
    {
      Gibbs(o);
    }
  }
}

//Equilibrazione Sistema
void Equilibrazione(int nstepeq)
{
  for(int i=1; i<=nstepeq; i++){
    Move(metro);
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  //int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = pow(u,2);
  walker[im] = m;
  walker[ix] = pow(m,2);
}


void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

   ofstream Ene, Heat, Mag, Chi;

    cout << "Block number " << iblk << endl;
//    cout << "Acceptance rate " << accepted/attempted << endl << endl;

//    cout << "---------------  en  -------------" << endl << endl;

    Ene.open("outputen.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << " " << setprecision(12) << iblk << " " <<  setprecision(12) << stima_u << " " << setprecision(12) << glob_av[iu]/(double)iblk << " " << setprecision(12) << err_u << endl;
    Ene.close();

    if(iblk == nblk){
      Ene.open("energie.0", ios::app);
      Ene << " " << temp << " " << setprecision(12) << glob_av[iu]/(double)iblk << " " << setprecision(12) << err_u << endl;
      Ene.close();
    }
//    cout << "---------------  cv  -------------" << endl << endl;

    Heat.open("outputcv.0", ios::app);
    stima_u2 = blk_av[ic]/blk_norm/(double)nspin;
    stima_c = pow(beta,2) * (stima_u2*(double)nspin - pow(stima_u,2) * (double)nspin*nspin)/(double)nspin;
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic], glob_av2[ic], iblk);
    Heat << " " << iblk << " " <<  setprecision(12) << stima_c << " " << setprecision(12) << glob_av[ic]/(double)iblk << " " <<  setprecision(12) << err_c << endl;
    Heat.close();

    if(iblk == nblk){
      Heat.open("cv.0", ios::app);
      Heat << " " << temp << " " << setprecision(12) << glob_av[ic]/(double)iblk << " " << setprecision(12) << err_c << endl;
      Heat.close();
    }
//    cout << "-------------  chi  ---------------" << endl << endl;

    Chi.open("outputchi.0", ios::app);
    stima_x = beta * blk_av[ix]/blk_norm/(double)nspin;
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix], glob_av2[ix], iblk);
    Chi << " " << iblk << " " <<  setprecision(12) << stima_x << " " << setprecision(12) << glob_av[ix]/(double)iblk << " " <<  setprecision(12) << err_x << endl;
    Chi.close();

    if(iblk == nblk){
      Chi.open("chi.0", ios::app);
      Chi << " " << temp << " " << setprecision(12) << glob_av[ix]/(double)iblk << " " << setprecision(12) << err_x << endl;
      Chi.close();
    }
//    cout << "------------  mag  ----------------" << endl << endl;

    Mag.open("outputmag.0", ios::app);
    stima_m =blk_av[im]/blk_norm/(double)nspin;
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im], glob_av2[im], iblk);
    Mag << " " << iblk << " " <<  setprecision(12) << stima_m << " " << setprecision(12) << glob_av[im]/(double)iblk << " " <<  setprecision(12) << err_m << endl;
    Mag.close();

    if(iblk == nblk){
      Mag.open("mag.0", ios::app);
      Mag << " " << temp << " " << setprecision(12) << glob_av[im]/(double)iblk << " " << setprecision(12) << err_m << endl;
      Mag.close();
    }
}

void Magnetizzazzione2(int iblk)
{
//  cout << "---------------  Mag2  -------------" << endl << endl;
  ofstream Mag2;
  Mag2.open("outputmag2.0", ios::app);
  stima_m = blk_av[im]/blk_norm/(double)nspin;
  glob_av[im] += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im], glob_av2[im], iblk);
  Mag2 << " " << iblk << setprecision(12) << stima_m << " " << setprecision(12) << glob_av[im]/(double)iblk << " " << setprecision(12) << err_m << endl;
  Mag2.close();

  if(iblk == nblk){
    Mag2.open("mag2.0", ios::app);
    Mag2 << " " << temp << " " << setprecision(12) << glob_av[im]/(double)iblk << " " << setprecision(12) << err_m << endl;
    Mag2.close();
  }
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void OldConfig(void)
{
  ofstream WriteConf;

  cout << " salva la configurazione in old.0" << endl << endl;
  WriteConf.open("old.0");
  for(int i=0; i<nspin; i++){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Metropolis(int o)
{
  double E = 2 * J * s[o] * (s[Pbc(o-1)] + s[Pbc(o+1)]) + h * 2. * s[o];
  if(E < 0) {
    s[o] = s[o] * (-1);
    accepted++;
  }
  else{
    double a = rnd.Rannyu();
    double prob = exp(-beta*E);
    if(a <= prob){
      s[o] = s[o] * (-1);
      accepted++;
    }
  }
  attempted++;
}

void Gibbs(int o)
{
  double E = 2*J*(s[Pbc(o-1)] + s[Pbc(o+1)]) + 2*h;
  double prob = 1/(1+exp(-beta * E));
  double a = rnd.Rannyu();

  if(a <= prob){
    s[o] = 1;
  }
  else{
    s[o] = -1;
  }
}
