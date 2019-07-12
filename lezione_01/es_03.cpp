#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

int main(){

// ***********************************************

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

double d=5.;					//spessore della piastrella circa 5cm
double L=3.5;					//lunghezza dell'ago circa 3.5cm
int M=100000;				//numero totale di numeri casuali usati
int N=100; 					//numero di blocchi
int m=int(M/N);				//numero di lanci per ogni blocco ------> E' FISSO!!!!!!
double pi_medio=0;
double pi2_medio=0;		//serve per dev std della media
double errore=0;

ofstream output("Data/Buffon.dat");

for(int i=0; i<N; i++){												//divido negli N blocchi
	double r,x,y,R;															//r unif tra [0,d/2], angolo theta ricavato da x e y 
																							//NB: sono x e y ad essere uniformi, non theta!!!		
	int N_hit=0;
	double pi;
	for(int n=0; n<m; n++){											//estraggo m numeri, se ago interseca aumento il contatore N_hit
		r=rnd.Uniform(0,d/2);
		if(r<L/2.){																//se punto medio dell'ago è più distante di L/2 dal bordo non lo interseca...
			do{
				x=rnd.Rannyu();												//x e y unif. all'interno di un quarto di cerchio unitario (altrimenti reject)
				y=rnd.Rannyu();
				R=sqrt(x*x + y*y);
			}while(R>=1);
			//theta=acos(x/R); ----> a noi interessa il cos(theta)
		
			if((x/R)>2*r/L)	N_hit++;								//-----> ago interseca la linea se: cos(theta)<L/(2*r) !!!!!
		}
	}
	pi=2*L*m/(double)(d*N_hit);													//valore di pi risultante da un blocco
	//cout << endl << pi;
	
	pi_medio=((pi_medio*i)+pi)/(double)(i+1);						//pi medio all'avanzare dei blocchi
	pi2_medio=((pi2_medio*i)+pi*pi)/(double)(i+1);
	
	if(i==0)	errore=0;
	else	errore=sqrt((pi2_medio-pow(pi_medio,2))/(double)i);			//dev std della media (all'avanzare dei blocchi)
	
	output << i+1 << " " << pi_medio << " " << errore << endl;
}
	
output.close();
return 0;
}
