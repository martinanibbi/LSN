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


int M=100000;			//numero totale di numeri casuali usati
int N=100; 				//numero di blocchi
int L=int(M/N);			//numero di lanci per ogni blocco ------> E' FISSO!!!!!!

ofstream output("Data/variance.dat");

double r_medio=0;
double r2_medio=0;
double errore=0;

for(int i=0; i<N; i++){																			//divido negli N blocchi
	double x=0, x2=0, y=0;																		//sommo numeri per media (per il singolo blocco)								
	
	for(int j=0; j<L; j++){																		//L lanci per ogni blocco, sommo i numeri generati
		y= rnd.Rannyu();																				//--->devo calcolare valore medio e dev.std di (r-1/2)^2
		y= (y-0.5)*(y-0.5);		
		x+= y;
	}
	
	double ave= x/(double)L;																						//media per singolo blocco
	double ave2= ave*ave;																					//media del quadrato per singolo blocco
	
	r_medio = ((r_medio*i)+ave)/(double)(i+1);									//media totale, all'avanzare dei blocchi
	r2_medio = ((r2_medio*i)+ave2)/(double)(i+1);

	if(i==0)	errore=0;																				//se un blocco solo non posso determinare l'incertezza, posta a 0!
		
	else	errore= sqrt((r2_medio-r_medio*r_medio)/(double)i);

	output << i+1 << " " << r_medio << " " << errore << endl;				//salva il numero dei blocchi, il valore medio e l'errore su un file
}


output.close();

return 0;
}

