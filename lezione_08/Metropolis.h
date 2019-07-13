#include "random.h"

using namespace std;

#ifndef _Metropolis_
#define _Metropolis_

double mu, sigma, delta;
double Emin = -0.46046587969271563;
		
double psi_trial(double x);		
	//densita' di probabilita' nello stato di prova psi_trial 
						
double E_loc(double x);
												
double step_psi_trial(double x, int*accept, Random* rnd);		
	//se cambia posizione (*accept)++
	
#endif
