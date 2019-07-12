#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "class_walk.h"

using namespace std;

Walk :: Walk(){
	r=0;
	n_steps=0;
}

Walk :: ~Walk(){}

double Walk :: distance(){
	return r;
}

// **********************

Walk_discrete :: Walk_discrete():Walk(){
	for(int i=0; i<3; i++)	posizione[i]=0;
}

Walk_discrete :: ~Walk_discrete(){
}

void Walk_discrete :: step(Random* rnd){
	double num=rnd->Rannyu();
	int indice;
	for(indice=1; indice<=6; indice++){					//divido intervallo [0,1] per assegnare la direzione e il verso:
																							//{1,3,5} indietro, {2,4,6} avanti, {1,2} x, {3,4} y, {5,6} z							 
		if(((double)indice/6.) > num)	break;				//indice---> numero intervallo 	=> 	indice/6 ------> bordo dx dell'intervallo
	}																						//mi fermo quando l'indice supera il numero
	int back_forward=1;
	if(indice%2!=0){
		back_forward=-1;																	//se indice numero dispari torno indietro, altrimenti vado avanti
		indice++;																					//porto {1,3,5} in {2,4,6}
	}
	indice = indice/2-1;																//porto {2,4,6} in {0,1,2}
	posizione[indice] += back_forward;
	r_value();																					//NB: alla fine di ogni step devo modificare r!!!
}

void Walk_discrete :: r_value(){
	double r2 = 0.;
	for(int i=0; i<3; i++)	r2 += (posizione[i])*(posizione[i]);
	r=sqrt(r2);
}

// **********************

Walk_continous :: Walk_continous():Walk(){
	for(int i=0; i<3; i++)	posizione[i]=0.;
}

Walk_continous :: ~Walk_continous(){
}

void Walk_continous :: step(Random* rnd){
	double phi=rnd->Uniform(0,2*M_PI);
	double T=rnd->Uniform(0,1);
	double theta=acos(1.-2.*T);	
	
	posizione[0]=posizione[0]+sin(theta)*cos(phi);
	posizione[1]=posizione[1]+sin(theta)*sin(phi);
	posizione[2]=posizione[2]+cos(theta);
	r_value();
}


void Walk_continous :: r_value(){
	double r2=0;
	for(int i=0; i<3; i++)	r2 += (posizione[i])*(posizione[i]);
	r=sqrt(r2);
}
