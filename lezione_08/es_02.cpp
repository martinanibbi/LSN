#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Metropolis.h"
#include "random.h"

using namespace std;

int main(){
//	leggo da Data/final_energy.dat i parametri: mu, sigma, energia
//	obiettivo: trovare mu e sigma che minimizzano l'energia
//	attenzione: energia non può scendere sotto Emin = -0.46046587969271563

	ifstream file_input;		//attenzione: se final_energy non viene svuotato ad ogni esecuzione di execute.py il codice potrebbe diventare 	
	file_input.open("Data/final_energy.dat");		//inefficente 
	//mu = 0.;
	//sigma = 0.;
	double energy_GS = 0.;		//energia del GS è sicuramente < 0
	int count = 0;	
	
	while(!file_input.eof() || count <= 10000){	//leggo tutto il file (count solo di controllo, dovrebbero esserci circa 100 righe)
		double mu_trial, sigma_trial, energy_trial;
		file_input >> mu_trial;
		file_input >> sigma_trial;
		file_input >> energy_trial;
		if(energy_trial <= energy_GS){
			energy_GS = energy_trial;
			sigma = sigma_trial;
			mu = mu_trial;
		}
		count ++;
	}
	
	ofstream final_GS;
	final_GS.open("Data/final_GS.dat");
	final_GS << mu << endl;
	final_GS << sigma << endl;
	final_GS << energy_GS << endl;
	if(energy_GS < Emin){
		cerr << endl << "ERRORE: l'energia del GS non può scendere sotto " << Emin << endl;
		final_GS << endl << endl << endl << "# ERRORE: l'energia del GS non può scendere sotto " << Emin << endl;
	}
	
file_input.close();
final_GS.close();

// ********************** campiono |psi_trial(x)|^2 **********************

// ************************
	 Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

// ********************************************


	//leggo posizione iniziale e delta da file
	double x;
	ifstream data_input;
	data_input.open("Data/input.dat");	
	data_input >> x;
	data_input >> delta;
	data_input.close();
	
	//double x = 1.;
	//delta = 0.2;
	//cout << endl << "mu: " << mu << endl << "sigma: " << sigma << endl;
	
	ofstream file_psi;
	file_psi.open("Data/psi_sampling.dat");
	int M = 100000;
	
	int accept=0;										//numero di valori accettati per ogni blocco (per controllo)
	for(int i=0; i < M; i++){				//eseguo il blocco
		x = step_psi_trial(x,&accept,&rnd);
		file_psi << x << endl;		
	}
	//cout << "Accettazione: " << (double)(100.*accept)/(double)M << "%" << endl;

file_psi.close();
return 0;
}

// *********************************************************************************** //


double psi_trial(double x){
	double e1 = exp(-(x-mu)*(x-mu)/(2.*sigma*sigma));
	double e2 = exp(-(x+mu)*(x+mu)/(2.*sigma*sigma));	
	double amp = e1 + e2;							
	return amp*amp;													
}

double step_psi_trial(double x, int*accept, Random* rnd){
	double new_x = x;
	new_x += rnd->Rannyu(-delta, delta);
	
	double A = min(1., psi_trial(new_x)/psi_trial(x));
	double m = rnd->Rannyu();
	
	if(m<A){
		x=new_x;
		(*accept)++;
	}	
	return x;
}

/*
double E_loc(double x){
	double s2 = sigma*sigma;
	double kin = (x*x + mu*mu)/s2 - 1. - (2.*x*mu/s2)*tanh(x*mu/s2);
	kin /= -2.*s2;		//ATTENZIONE: termine cinetico è negativo per (-1)^2 e moltiplicato per h^2/2m (ma h e m == 1)
	double pot = (x*x - 2.5)*x*x;
	return kin + pot;
}
*/		//non serve

