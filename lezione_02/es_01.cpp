#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

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

//*************** (1) *****************
//distribuzione uniforme, metodo della media

int M=100000;										//numero totale di numeri estratti
int N=100;											//numero totale di blocchi
int L=M/N;											//numero di lanci per ogni blocco
double I_medio=0;
double I2_medio=0;
double errore=0;

ofstream out_unif("Data/integrale_uniforme.dat");

for(int n=0; n<N; n++){					//aumento il numero di blocchi
	double f=0;

	for(int l=0; l<L; l++){				//estraggo x uniformemente L volte e mi calcolo la media delle f(x) che ottengo ----> I=f_medio*intervallo
		double x=rnd.Rannyu();
		f+=integrand(x);
	}
	f /= (double)L;
	I_medio = ((I_medio*n) + f)/(double)(n+1);	//faccio UN calcolo dell'integrale per ogni blocco e ne faccio la media all'aumentare dei blocchi
	I2_medio=((I2_medio*n) + f*f)/(n+1);
	
	if(n!=0)	errore=sqrt((I2_medio-I_medio*I_medio)/n);					//errore=dev std della media, al primo passaggio errore=0	
	
	out_unif << n+1 << " " << I_medio << " " << errore << endl;
}

out_unif.close();

//*************** (2) *****************
//distribuzione: p(x)=2(1-x) (gia' normalizzata), metodo della media
//con metodo funzione inversa: x=1-sqrt(1-z) con z uniforme ------> scelgo sol (-) per avere p(x) tra [0,2] => x tra [0,1]

I_medio=0;
I2_medio=0;
errore=0;

ofstream out_samp("Data/integrale_sampling.dat");

for(int n=0; n<N; n++){					//aumento il numero di blocchi
	double f=0;

	for(int l=0; l<L; l++){				//estraggo x uniformemente L volte e mi calcolo la media delle f(x) che ottengo ----> I=f_medio*intervallo
		double x=rnd.Rannyu();
		double y=1.-sqrt(1-x);			//distribuzione ottenuta con funzione inversa
		f+=new_integrand(y);
	}
	f /= (double)L;
	I_medio = ((I_medio*n)+f)/(double)(n+1);	//faccio UN calcolo dell'integrale per ogni blocco e ne faccio la media all'aumentare dei blocchi
	I2_medio=((I2_medio*n)+ f*f)/(double)(n+1);
	
	if(n!=0)	errore=sqrt((I2_medio-I_medio*I_medio)/(double)n);					//errore=dev std della media, al primo passaggio errore=0	
	
	out_samp << n+1 << " " << I_medio << " " << errore << endl;
}

out_samp.close();

return 0;
}

double integrand(double x){
	return M_PI*cos(M_PI*x/2.)/2.;
}

double distribution(double x){
	return 2.*(1-x);
}

double new_integrand(double x){									//non dovrebbe esserci il rischio di avere /0 x' x tra [0,1]
	return integrand(x)/distribution(x);
}
