#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
#include "Travelling_Salesman.h"

using namespace std;

int main(){

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
	
	Sim_Ann SA(&rnd,10);
	SA.Print_City();
	
	for(double i=10; i>1; i-=0.05){
		SA.Set_Temp(i);
		for(int j=0; j<1000; j++){
			SA.Step();
		}
		SA.Print_Path();	//stampo miglior elemento dopo 10000 step 
	}
	for(double i=1; i>0; i-=0.005){
	SA.Set_Temp(i);
		for(int j=0; j<10000; j++){
			SA.Step();
		}
		SA.Print_Path();
	}
	for(double i=0.0049; i>0.0001; i-=0.0001){
	SA.Set_Temp(i);
		for(int j=0; j<10000; j++){
			SA.Step();
		}
		SA.Print_Path();
	}
	SA.Print_Conf();	
	
return 0;
}

