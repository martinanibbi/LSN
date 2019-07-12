#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "random.h"
#include "Metropolis.h"

using namespace std;

double state_1s(double X, double Y, double Z){
	double r_a0 = R(X,Y,Z);
	double amp = exp(-r_a0);								//ampiezza di probabilita' dello stato fondamentale
	return amp*amp;													//NB: manca un fattore (1/pi), ma tanto nel calcolo di A si elide...
}

double state_2p(double X, double Y, double Z){
	double r_a0 = R(X,Y,Z);
	double amp = r_a0*exp(-r_a0/2.)*cos_z(X,Y,Z);
	return amp*amp;
}

double state_3d(double X, double Y, double Z){
	double r_a0 = R(X,Y,Z);
	double amp = r_a0*r_a0*exp(-r_a0/3.)*(3*pow(cos_z(X,Y,Z),2)-1);
	return amp*amp;
}

double state_4f(double X, double Y, double Z){
	double r_a0 = R(X,Y,Z);
	double amp = r_a0*r_a0*r_a0*exp(-r_a0/4.)*(5*pow(cos_z(X,Y,Z),3)-3*cos_z(X,Y,Z));
	return amp*amp;
}

position Step(position pos, int*accept){
	if(distribution == "Uniform"){
		if(state == "1s")	return Step_1s_unif(pos, accept);
		if(state == "2p")	return Step_2p_unif(pos, accept);
		if(state == "3d")	return Step_3d_unif(pos, accept);
		if(state == "4f")	return Step_4f_unif(pos, accept);
	}
	if(distribution == "Gauss"){
		if(state == "1s")	return Step_1s_gaus(pos, accept);
		if(state == "2p")	return Step_2p_gaus(pos, accept);
		//if(state == "3d")	return Step_3d_gaus(pos, accept);
		//if(state == "4f")	return Step_4f_gaus(pos, accept);	 //NB: per Gauss implementati solo stati 1s e 2p...
	}
}


position Step_1s_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd->Uniform(-delta, delta);
	new_pos.y += rnd->Uniform(-delta, delta);
	new_pos.z += rnd->Uniform(-delta, delta);
	
	double A = min(1., state_1s(new_pos.x,new_pos.y,new_pos.z)/state_1s(pos.x,pos.y,pos.z));
	double m = rnd->Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_1s_gaus(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd->Gauss(0, delta/2.);
	new_pos.y += rnd->Gauss(0, delta/2.);
	new_pos.z += rnd->Gauss(0, delta/2.);
	
	double A = min(1., state_1s(new_pos.x,new_pos.y,new_pos.z)/state_1s(pos.x,pos.y,pos.z));
	double m = rnd->Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_2p_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd->Uniform(-delta, delta);
	new_pos.y += rnd->Uniform(-delta, delta);
	new_pos.z += rnd->Uniform(-delta, delta);
	
	double A = min(1., state_2p(new_pos.x,new_pos.y,new_pos.z)/state_2p(pos.x,pos.y,pos.z));
	double m = rnd->Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_2p_gaus(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd->Gauss(0, delta/2.);
	new_pos.y += rnd->Gauss(0, delta/2.);
	new_pos.z += rnd->Gauss(0, delta/2.);
	
	double A = min(1., state_2p(new_pos.x,new_pos.y,new_pos.z)/state_2p(pos.x,pos.y,pos.z));
	double m = rnd->Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_3d_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd->Uniform(-delta, delta);
	new_pos.y += rnd->Uniform(-delta, delta);
	new_pos.z += rnd->Uniform(-delta, delta);
	
	double A = min(1., state_3d(new_pos.x,new_pos.y,new_pos.z)/state_3d(pos.x,pos.y,pos.z));
	double m = rnd->Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_4f_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd->Uniform(-delta, delta);
	new_pos.y += rnd->Uniform(-delta, delta);
	new_pos.z += rnd->Uniform(-delta, delta);
	
	double A = min(1., state_4f(new_pos.x,new_pos.y,new_pos.z)/state_4f(pos.x,pos.y,pos.z));
	double m = rnd->Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

double R(double x, double y, double z){
	return sqrt( x*x + y*y + z*z );
}

double cos_z(double x, double y, double z){
	return z/R(x,y,z);
}
