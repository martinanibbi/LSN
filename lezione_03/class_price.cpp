#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "class_price.h"

using namespace std;

Price :: Price(double s0, double t, double k, double R, double Sigma){
	S0=s0;
	T=t;
	K=k;
	r=R;
	sigma=Sigma;
}


Price :: ~Price(){
}

double Price :: Call_price(Random* rnd){
	double S=asset_price(rnd);
	if(S<=K)	return 0;	
	else	return (S-K)*exp(-r*T);												//NB: se S-K e' negativo, profitto e' nullo.
}

double Price :: Put_price(Random*rnd){
	double S=asset_price(rnd);
	if(S>=K)	return 0;	
	else	return (K-S)*exp(-r*T);												//NB: se K-S e' negativo, profitto e' nullo.
}


// *************** CONTINUO ***************

Price_direct :: Price_direct(double s0, double t, double k, double R, double Sigma) : Price(s0, t, k, R, Sigma){
}

Price_direct :: ~Price_direct(){
}

double Price_direct :: asset_price(Random* rnd){
	double W=rnd->Gauss(0,sqrt(T));											//attenzione: sigma^2=T !!!!s
	return S0*exp(T*(r-pow(sigma,2)/2.)+sigma*W);
}


// ************** DISCRETO ***************


Price_discrete :: Price_discrete(double s0, double t, double k, double R, double Sigma) : Price(s0, t, k, R, Sigma){
}

Price_discrete :: ~Price_discrete(){
}

double Price_discrete :: asset_price(Random* rnd){
	double S=S0;
	for(double t=0.01; t<=1.; t+=0.01){										//divido intervallo di tempo [0,1] in 100
		double W=rnd->Gauss(0.,1.);													//ogni volta genero un numero da Gaussiana a sigma fissata=1
		S=S*exp(0.01*(r-pow(sigma,2)/2.)+sigma*W*0.1);
	}
	return S;
}


