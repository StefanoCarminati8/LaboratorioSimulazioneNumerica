/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "random.h"

//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie, itot;
double stima_pot, stima_kin, stima_etot, stima_temp;

Random rnd;
int seed[4];

int nblocchi;

// averages
double acc,att;
double blk_av[m_props], blk_norm, accepted,attempted;
double glob_av[m_props], glob_av2[m_props];
double err_pot,err_temp,err_kin,err_etot, temperatura;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part], xold2[m_part],yold2[m_part],zold2[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double v, t;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut, scenario;

// simulation
int nstep, iprint;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Termalizza(void);
void Velocita(void);
void VelocitaQ(void);
void OldConfig();
void SuperOldConfig();
double Error(double,double,int);
void Reset(int);
void Accumula(void);
void Medie(int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
