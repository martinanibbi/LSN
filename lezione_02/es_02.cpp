#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "class_walk.h"

using namespace std;

double integrand(double x);
double distribution(double x);
double new_integrand(double x);

int main(){

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

// ***********************************************

int N_steps=100;
int N_blocchi=100;
int M=10000;													//numero punti totali estratti(per ogni step)
int L=M/N_blocchi;										//numero punti per ogni blocco
double r2_medio=0;										//devo eseguire la media del modulo^2 della distanza (attenzione, x calcolo varianza ci vorra'
double r4_medio=0;
double errore=0;

// ******************* discreto ****************

ofstream out_d("Data/random_discreto.dat");
ofstream prova("Data/prova.dat");
//ofstream scatter("Data/scatter.dat");

for(int i=1; i<=N_steps; i++){									//divido per numero di step, eseguo M simulazioni per ogni step
	
	for(int n=0; n<N_blocchi; n++){								//eseguo l'esperimento M volte, diviso in blocchi
		
		double r2_medio_blocco=0.;
		double r4_medio_blocco=0.;
		
		for(int l=0; l<L; l++){													//eseguo il blocco
			Walk* drunk = new Walk_discrete();
			for(int j=0; j<i; j++){
				drunk->step(&rnd);		//arrivo al passo i-esimo			NB: devo passargli rnd by reference!!
			}		
			r2_medio_blocco += (drunk->distance())*(drunk->distance());
			//if( n==(N_blocchi-1) && i==N_steps )	scatter << drunk->GetX() << " " << drunk->GetY() << " " << drunk->GetZ() << endl;
			delete drunk;
		}
		
		r2_medio_blocco /= (double)L;
		r4_medio_blocco = r2_medio_blocco*r2_medio_blocco;
		
		r2_medio=((r2_medio*n)+r2_medio_blocco)/(double)(n+1.);				//aggiorno la media ad ogni blocco
		r4_medio=((r4_medio*n)+r4_medio_blocco)/(double)(n+1.);
		if(i == N_steps)	prova << n+1 << " " << r2_medio << " " << sqrt((r4_medio-r2_medio*r2_medio)/n) << endl;
	}
	
	errore=(sqrt((r4_medio-r2_medio*r2_medio)/(N_blocchi-1.)))/(2.*sqrt(r2_medio));			//NB: ERRORE SULLA RADICE DI <r^2>!!!!!!
	out_d << i << " " << sqrt(r2_medio) << " " << errore << endl;
}

out_d.close();
prova.close();
//scatter.close();

// ******************* continuo *****************

ofstream out_c("Data/random_continuo.dat");

for(int i=1; i<=N_steps; i++){									//divido per numero di step, eseguo M simulazioni per ogni step
	
	for(int n=0; n<N_blocchi; n++){								//eseguo l'esperimento M volte, diviso in blocchi
		
		double r2_medio_blocco=0;
		double r4_medio_blocco=0;
				
		for(int l=0; l<L; l++){													//eseguo il blocco
			Walk* drunk = new Walk_continous();
			for(int j=0; j<i; j++)	drunk->step(&rnd);		//arrivo al passo i-esimo			NB: devo passargli rnd by reference!!
			r2_medio_blocco += (drunk->distance())*(drunk->distance());
			delete drunk;
		}
		
		r2_medio_blocco /= (double)L;
		r4_medio_blocco = r2_medio_blocco*r2_medio_blocco;
		
		r2_medio = ((r2_medio*n)+r2_medio_blocco)/(double)(n+1.);				//aggiorno la media ad ogni esperimento
		r4_medio = ((r4_medio*n)+r4_medio_blocco)/(double)(n+1.);
	}
	errore = sqrt((r4_medio-r2_medio*r2_medio)/(N_blocchi-1.))/(2.*sqrt(r2_medio));			//NB: ERRORE SULLA RADICE DI <r^2>!!!!!!
	out_c << i << " " << sqrt(r2_medio) << " " << errore << endl;
}

out_c.close();

return 0;
}


