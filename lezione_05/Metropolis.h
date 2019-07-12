#include "random.h"

using namespace std;

#ifndef _Metropolis_
#define _Metropolis_

//letti da input.dat
double x_0, y_0, z_0, delta;		
string state; 		//1s, 2p, 3d, 4f
string distribution; 	//Gauss, Uniform
string folder;	//Data/'stato'/'distribuzione'/
bool start;	//se false eseguo solo la fase di equilibrazione (serve per determinare delta)
Random rnd;

struct position{
	double x;
	double y;
	double z;
};
			
double state_1s(double X, double Y, double Z);		
	//densita' di probabilita' nello stato (1,0,0) 
	
double state_2p(double X, double Y, double Z);
	//densità di probabilità nello stato (2,1,m) con m = -1, 0, 1
	
double state_3d(double X, double Y, double Z);
	//densità di probabilità nello stato (3,2,m) con m = -2, -1, 0, 1, 2
	
double state_4f(double X, double Y, double Z);
	//densità di probabilità nello stato (4,3,m) con m = -3, -2, -1, 0, 1, 2, 3
						
position Step(position pos, int*accept);	
	//rimanda a uno Step specifico in base allo stato e alla distribuzione letti da input.dat											
position Step_1s_unif(position pos, int*accept);		
position Step_2p_unif(position pos, int*accept);	
position Step_3d_unif(position pos, int*accept);			
position Step_4f_unif(position pos, int*accept);			
position Step_1s_gaus(position pos, int*accept);		
position Step_2p_gaus(position pos, int*accept);		
	//se cambia posizione (*accept)++
	
double R(double X, double Y, double Z);			
	//distanza dall'origine
	
double cos_z(double X, double Y, double Z);
	//cos(theta) con momento angolare diretto lungo l'asse z 

#endif
