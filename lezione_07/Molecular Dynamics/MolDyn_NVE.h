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

using namespace std;

Random rnd;

//pigreco
const double pi=3.1415927;

//preliminary input
bool real;																	//true ---> fluido reale (A/K);		//false ---> unità di Lennard-Jones
bool frames;
string phase, element, folder, file_path;		//folder distingue solo tra s, l, g		//filepath distingue ulteriormente tra A, K, LJU

//parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblk;
double delta;

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,err_pot,err_press,err_gdir;

//Real fluid
double m, sigma, e_kb, epsilon;							 
double kb = 1.38064852e-23;									//per conversione in unità di misura del SI mi servono costante di Boltzmann,
double amu = 1.66054e-24;										//unità di massa atomica,
double nm = 1e-9;														//nanometri
																						//rho, energia, vol etc. reali sovrascritti su quelli calcolati in unità di LJ
//functions
void Random_start(void);
void Input(void);
	void Input_phase(string);
	void Input_conf(int);
		void Input_conf_FCC(void);
		void Input_conf_old(void);
	void Input_element(string);
	
void Move(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
