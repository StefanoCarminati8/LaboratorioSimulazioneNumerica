/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include <vector>
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

int main(){
  Input();             //Inizialization

  cout << "temperatura attuale: " << temp << endl;
  ofstream Term;
  Term.open("output_restart.out");

  int steps = 20000;
  cout << "termalizzazione del sistema" << endl;
  for(int i=0; i<steps;i++){
    Move();

    if(i%10 == 0){
      Termalizza();
      Term << stima_temp << endl;
    }
  }
  Term.close();

  cout << "temperatura post termalizzazione: " << stima_temp << endl;
  cout << "                                            " << endl;

  cout << "inizio simulazione" << endl;

  int nconf = 1;
  for(int iblk = 1; iblk<nblocchi; iblk++){
    cout << "blocco numero: " << iblk << endl;
    Reset(iblk);
    for(int istep=1; istep <= nstep; ++istep){
       Move();           //Move particles with Verlet algorithm
       if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
       if(istep%10 == 0){
          Measure();     //Properties measurement
  //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
          Accumula();
          nconf += 1;
       }
    }
    Medie(iblk);
  }
  OldConfig();
  SuperOldConfig();
  ConfFinal();         //Write final configuration to restart
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
//  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //Preparo generatore numeri casuali
  int p1, p2=0;
  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  ifstream seeds("seed.in");
  seeds >> seed[0] >> seed[1]  >> seed[2]  >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  seeds.close();

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblocchi;
  ReadInput >> scenario;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //energy
  it = 3; //Temperature
  itot = 4; //energia totale
  n_props = 4; //Number of observables

//Read initial configuration
  if(scenario == 0){ //primo step
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
  }

  if(scenario == 1){ //secondo step
    cout << "Read old configuration from file old.0 " << endl << endl;
    ReadConf.open("old.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      x[i] = xold[i] * box;
      y[i] = yold[i] * box;
      z[i] = zold[i] * box;
    }
    ReadConf.close();
  }

  if(scenario == 2){ //secondo step
    cout << "Read extra old configuration from file old.final " << endl << endl;
    ReadConf.open("old.final");
    for (int i=0; i<npart; ++i){
      ReadConf >> xold2[i] >> yold2[i] >> zold2[i];
      x[i] = xold2[i] * box;
      y[i] = yold2[i] * box;
      z[i] = zold2[i] * box;
    }
    ReadConf.close();
  }
//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i) {
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

void Termalizza(void){
  double sumv[3] = {0.0, 0.0, 0.0};
  for(int i=0; i<npart; i++){
    vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);
    vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
    vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for(int dim=0; dim<3; dim++){
    sumv[dim] /= (double)npart;
  }
  double sumv2 = 0.0;
  double fs=0.0;
  for(int i=0; i<npart; i++){
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];
    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;

  VelocitaQ();
  temp = (2/3) * t / (double)npart;
  stima_temp = temp;

  fs = sqrt(3*temp / sumv2);
  for(int i=0; i<npart; i++){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;
  }
  return;
}

void Velocita(void){
  for(int i=0; i<npart; i++){
    vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);
    vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
    vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
  }
  return;
}

void VelocitaQ(void){
  double t=0.;
  for(int i=0; i<npart; i++){
    t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xold2[i] = xold[i];
    yold2[i] = yold[i];
    zold2[i] = zold[i];


    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(){ //Properties measurement
//  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void OldConfig(void){
  ofstream oldconfig;

  cout << "Salva vecchia configurazione in old.0 " << endl << endl;
  oldconfig.open("old.0");
  for(int i=0; i<npart; i++){
    oldconfig << xold[i]/box << "   " << yold[i]/box << "   " << zold[i]/box << endl;
  }
  oldconfig.close();
}

void SuperOldConfig(void){
  ofstream superold;
  cout << "Salva configurazione super vecchia in old.final " << endl << endl;
  superold.open("old.final");

  for(int i=0; i<npart; i++){
    superold << xold2[i]/box << "   " << yold2[i]/box << "   " << zold2[i]/box << endl;
  }
  superold.close();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk){
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk, 2))/(double)iblk);
}

void Reset(int iblk){
  if(iblk == 1){
    for(int i=0; i<n_props; i++){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }
  for(int i=0; i<n_props; i++){
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumula(void){
  for(int i=0; i<n_props; i++){
    blk_av[i] = blk_av[i] * delta;
  }
}

void Medie(int iblk){
  ofstream epot, temp, ekin, etot;
  epot.open("media_epot.out");
  temp.open("media_temp.out");
  ekin.open("media_ekin.out");
  etot.open("media_etot.out");

  //energia potenziale
  stima_pot = blk_av[iv]/blk_norm/(double)npart;
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot = Error(glob_av[iv], glob_av2[iv], iblk);

  //temperatura
  temperatura = blk_av[it]/blk_norm;
  glob_av[it] += temperatura;
  glob_av2[it] += temperatura*temperatura;
  err_temp = Error(glob_av[it], glob_av2[it], iblk);

  //energia cinetica
  stima_kin = blk_av[ik]/blk_norm;
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin = Error(glob_av[ik], glob_av2[ik], iblk);

  //energia totale
  stima_etot = stima_kin + stima_pot;
  glob_av[itot] += stima_etot;
  glob_av2[itot] += stima_etot*stima_etot;
  err_etot = Error(glob_av[itot], glob_av2[itot], iblk);

  //valori per particella
  epot << " " << iblk << setprecision(12) << glob_av[iv]/(double)iblk << " " << setprecision(12) << err_pot << endl;
  temp << " " << iblk << setprecision(12) << glob_av[it]/(double)iblk << " " << setprecision(12) << err_temp << endl;
  ekin << " " << iblk << setprecision(12) << glob_av[ik]/(double)iblk << " " << setprecision(12) << err_kin << endl;
  etot << " " << iblk << setprecision(12) << glob_av[itot]/(double)iblk << " " << setprecision(12) << err_etot << endl;

  cout << endl << " _________calcolando...._________" << endl << endl;
  epot.close();
  temp.close();
  ekin.close();
  etot.close();
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
