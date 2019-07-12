#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "class_price.h"

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

double S0 = 100.;												//prezzo a t=0
double K = 100.;												//strike price
double T = 1.;													//tempo finale
double r = 0.1;													//tasso di interesse 
double sigma = 0.25;										//volatilita'

int M=100000;														//numero di simulazioni (?)
int N=100;															//numero di blocchi
int L=M/N;															//numero di simulazioni per blocco

// ******************* CONTINUO *********************

double Call_medio=0, Call2_medio=0;
double Put_medio=0, Put2_medio=0;
double errore_Call=0, errore_Put=0;

ofstream CallC_out("Data/Call_direct.dat");
ofstream PutC_out("Data/Put_direct.dat");

Price* p_c= new Price_direct(S0,T,K,r,sigma);

for(int n=0; n<N; n++){									//divido nei blocchi
	double Call=0, Call2=0;
	double Put=0, Put2=0;
	
	for(int l=0; l<L; l++){								//in ogni blocco faccio la media dei prezzi
		double x= p_c->Call_price(&rnd);
		double y= p_c->Put_price(&rnd);
		Call+= x;
		Put+= y;
	}
	
	Call = double(Call/L);									//media su un blocco
	Call2 = Call*Call;
	Put = double(Put/L);
	Put2 = Put*Put;
	
	Call_medio= ((Call_medio*n)+Call)/(n+1);
	Call2_medio= ((Call2_medio*n)+Call2)/(n+1);
	Put_medio= ((Put_medio*n)+Put)/(n+1);
	Put2_medio= ((Put2_medio*n)+Put2)/(n+1);
	
	if(n>0){
		errore_Call=sqrt( (Call2_medio - Call_medio*Call_medio)/double(n));
		errore_Put=sqrt( (Put2_medio - Put_medio*Put_medio)/double(n));
	}
	
	CallC_out << n+1 << " " << Call_medio << " " << errore_Call << endl;
	PutC_out << n+1 << " " << Put_medio << " " << errore_Put << endl;
}

CallC_out.close();
PutC_out.close();

// ******************** DISCRETO *******************


Call_medio=0;
Call2_medio=0;
Put_medio=0;
Put2_medio=0;
errore_Call=0;
errore_Put=0;

ofstream CallD_out("Data/Call_discrete.dat");
ofstream PutD_out("Data/Put_discrete.dat");

Price* p_d= new Price_discrete(S0,T,K,r,sigma);

for(int n=0; n<N; n++){									//divido nei blocchi
	double Call=0, Call2=0;
	double Put=0, Put2=0;
	
	for(int l=0; l<L; l++){								//in ogni blocco faccio la media dei prezzi
		double x= p_d->Call_price(&rnd);		
		double y= p_d->Put_price(&rnd);			
		Call+= x;
		Put+= y;
	}
	
	Call = double(Call/L);									//media su un blocco
	Call2 = Call*Call;
	Put = double(Put/L);
	Put2 = Put*Put;
	
	Call_medio= ((Call_medio*n)+Call)/(n+1);
	Call2_medio= ((Call2_medio*n)+Call2)/(n+1);
	Put_medio= ((Put_medio*n)+Put)/(n+1);
	Put2_medio= ((Put2_medio*n)+Put2)/(n+1);
	
	if(n>0){
		errore_Call=sqrt( (Call2_medio - Call_medio*Call_medio)/double(n));
		errore_Put=sqrt( (Put2_medio - Put_medio*Put_medio)/double(n));
	}
	
	CallD_out << n+1 << " " << Call_medio << " " << errore_Call << endl;
	PutD_out << n+1 << " " << Put_medio << " " << errore_Put << endl;
}

CallD_out.close();
PutD_out.close();

return 0;
}
