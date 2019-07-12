/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <string>

using namespace std;

#ifndef __Ising_
#define __Ising_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig;
double nbins;
double walker[m_props];					

/*
	In walker salvo valori "instantanei" (ma di tutta la configurazione) con Measure, in particolare:
	0 -> u   (da cui U = <u>)
	1 -> u^2 (da cui C ~ <u^2> -<u>^2)
	2 -> m^2 (da cui X ~ <m^2>)
	3 -> m   (da cui M = <m>)
*/

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x,stima_g;
double err[m_props];

//configuration
const int m_spin=50;
double s[m_spin];

// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, nblk;
string metro;
string filename;		//nome file dipende dalla temperatura e campo magnetico, file divisi in cartelle Metropolis/Gibbs e Energy/Pressure etc.

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move();
void ConfFinal(void);
void Measure(void);
double Boltzmann(int sm, int ip);								// sm ---> valore dello spin (+-1);		ip ---> posizione nell'array (chiamato come o...)
int Pbc(int);
double Error(double,double,int);

int Metro_Boltzmann(int o);											//prende in ingresso posizione della particella, restituisce 0 se scartato, 1 se accettato
void Gibbs_Boltzmann(int o);										//NB: valore sempre accettato...
void Config_Random();
void Config_Old();
void Equilibration();

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
