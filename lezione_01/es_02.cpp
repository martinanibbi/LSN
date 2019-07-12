#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

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

/*ofstream prova("prova");
for(int i=0; i<10000; i++){
	prova << rnd.Exponential(1) << endl;
}
prova.close();
*/																			//prova distribuzione esponenziale

int J=10000;														//numero di volte in cui eseguo la somma (riempio l'istogramma)

ofstream** output=new ofstream*[3];
	for(int i=0; i<3; i++)	output[i]=new ofstream[4];
		output[0][0].open("Data/unif_1.dat");
		output[0][1].open("Data/unif_2.dat");
		output[0][2].open("Data/unif_10.dat");
		output[0][3].open("Data/unif_100.dat");
		output[1][0].open("Data/exp_1.dat");
		output[1][1].open("Data/exp_2.dat");
		output[1][2].open("Data/exp_10.dat");
		output[1][3].open("Data/exp_100.dat");
		output[2][0].open("Data/lor_1.dat");
		output[2][1].open("Data/lor_2.dat");
		output[2][2].open("Data/lor_10.dat");
		output[2][3].open("Data/lor_100.dat");						//apro 3*4 files
		
int N[4];																
N[0]=1;
N[1]=2;
N[2]=10;
N[3]=100;

for(int i=0; i<4; i++){									//vario numero su cui sommo
	double* Sum=new double[3];						//tre tipi di distribuzioni
	
	for(int j=0; j<J; j++){								//eseguo l'esperimento J volte
		double Sum_0=0, Sum_1=0, Sum_2=0;
		for(int n=0; n<N[i]; n++){					//eseguo la somma
			Sum_0+=rnd.Rannyu();							// 0) UNIFORME
			Sum_1+=rnd.Exponential(1);				// 1) ESPONENZIALE	(lambda=1)
			Sum_2+=rnd.Lorentzian(1,0);				// 2) LORENTZIANA		(Lambda=1, mu=0)
		}
		Sum[0]=double(Sum_0/N[i]);
		Sum[1]=double(Sum_1/N[i]);
		Sum[2]=double(Sum_2/N[i]);
		
		for(int k=0; k<3; k++)	output[k][i] << Sum[k] << endl;		//stampo unif, exp, lor su 3 file diversi (a Ni fissato)
	}	
	delete[]Sum;
}

for(int i=0; i<3; i++){																//chiudo e cancello i vettori
	for(int j=0; j<4; j++)	output[i][j].close();
	delete[]output[i];
}
delete[]output;

return 0;
}

// PER JUPYTER: fare un fit con una Gaussiana per N=100 per unif e exp!!!!!!




