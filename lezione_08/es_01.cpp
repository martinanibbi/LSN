#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "Metropolis.h"

using namespace std;

int main(){
int wd = 10;

// ******************************
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

// ********************************
	
	//leggo posizione iniziale, delta, media e dev std da file
	double x;
	ifstream data_input;
	data_input.open("Data/input.dat");	
	data_input >> x;
	data_input >> delta;
	data_input >> mu;
	data_input >> sigma;

	cout << endl << endl << "Starting position: " << x;
	cout << endl << "Delta: " << delta;
	cout << endl << "Gaussian mean: " << mu;
	cout << endl << "Gaussian standard deviation: " << sigma << endl << endl;
	
// ********************************

	ofstream out_energy;
	ofstream final_energy;
	out_energy.open("Data/energy.dat");
	final_energy.open("Data/final_energy.dat", ios::app);	//salvo ultimo calcolo dell'energia con valori di mu e sigma
	
	int N_blocchi = 100;
	int M = 100000;
	int L = M/N_blocchi;																				//numero di passi per ogni blocco
	
	double average=0., average2=0., error=0.;
	for(int n=0; n<N_blocchi; n++){		//divido nei blocchi
		int accept=0;																							//numero di valori accettati per ogni blocco (per controllo)
		double block=0., block2=0.;
		for(int istep=0; istep < L; istep++){		//eseguo il blocco
			x = step_psi_trial(x, &accept, &rnd);
			block += E_loc(x);		
			//if(n==0)	cout << endl << x;
		}
		block = block/(double)L;
		block2 = block*block;
		
		//cout << "Blocco: " << n+1. << "/" << N_blocchi << endl;
		//cout << "Accettazione: " << (double)(100.*accept)/(double)L << "%" << endl;

		average = (double)(((average)*n + block)/(double)(n+1.));
  	average2 = (double)(((average2)*n + block2)/(double)(n+1.));				//medie all'avanzare dei blocchi
  	if( n!=0)	error = sqrt( (average2 - average*average)/(double)n);
		
		out_energy << n+1. << " " << setw(wd) << average << " " << setw(wd) << error << endl;
	}
	final_energy << setw(wd) << mu << " " << setw(wd) << sigma << " " << setw(wd) << average << endl;	//non salvo errore
	
	out_energy.close();
	final_energy.close();

return 0;  
}

// **********************************************************************************************//

double psi_trial(double x){
	double e1 = exp(-(x-mu)*(x-mu)/(2.*sigma*sigma));
	double e2 = exp(-(x+mu)*(x+mu)/(2.*sigma*sigma));	
	double amp = e1 + e2;							
	return amp*amp;	//normalizzazione si elide nel Metropolis												
}

double step_psi_trial(double x, int*accept, Random* rnd){
	double step = delta*(1-2.*rnd->Rannyu());
	double new_x = x + step;
	double A = min(1., psi_trial(new_x)/psi_trial(x));
	double m = rnd->Rannyu();
	
	if(m<A){
		x=new_x;
		(*accept)++;
	}	
	return x;
}

double E_loc(double x){
	double s2 = sigma*sigma;
	double kin = (x*x + mu*mu)/s2 - 1. - (2.*x*mu/s2)*tanh(x*mu/s2);
	kin /= (-2.*s2);		//ATTENZIONE: termine cinetico è negativo per (-i)^2 e moltiplicato per h^2/2m (ma h e m == 1)
	double pot = (x*x - 2.5)*x*x;
	return kin + pot;
}


