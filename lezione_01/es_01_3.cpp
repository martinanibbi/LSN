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

int M=100;										//divido [0,1] in M intervalli e riempio l'istogramma
int N=10000;									//numero totale di numeri casuali usati
int J=100;										//numero di volte in cui ripeto l'esperimento (calcolo J chi quadro)
ofstream output("Data/chi_squares.dat");

double E= N/M;									//valore medio e incertezza!
int*n= new int[M];							//n[i]= # di numeri nell'i-esimo intervallo

for(int j=0; j<J; j++){
	
	for(int i=0; i<M; i++)	n[i]=0;					//svuoto il vettore ad ogni esperimento
	
	for(int i=0; i<N; i++){			//riempio l'istogramma
		double y=rnd.Rannyu();		//genero un numero tra 0 e 1
						
															//m=estremo destro dell'intervallo moltiplicato per M!!! (indica quale n[i] far aumentare di 1)
		for(int m=1; m<=M; m++){	//faccio scorrere gli intervalli, appena y supera l'estremo dx mi fermo
			if(M*y<double(m)){
				n[m-1]++;
				//cout << endl << n[m-1];
				break;
			}
		}
	}
	double x=0;
	for(int k=0; k<M; k++){			//calcolo il chi quadro con il vettore n riempito in precedenza
		x+=pow((double(n[k])-E),2);	
	}
	x=x/E;
	output << j+1 << " " << x << endl;
	//cout << endl << x;
}
output.close();
delete[]n;
return 0;
}
