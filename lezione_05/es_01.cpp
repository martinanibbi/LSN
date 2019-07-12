#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "Metropolis.h"

using namespace std;

int main(){

// ************************

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

// *************** READ INPUT ******************
	
	ifstream ReadInput("Data/input.dat");	//leggo posizione iniziale, delta e stato da file 
	ReadInput >> x_0 >> y_0 >> z_0;
	position pos = {x_0, y_0, z_0};	//posizione iniziale
	ReadInput >> delta;	//delta
	ReadInput >> state;	//1s, 2p, 3d, 4f
	ReadInput >> distribution;	//Uniform, Gauss
	string answer;
	ReadInput >> answer;
	if(answer == "y") start = true;
	else{
		start = false;
		cout << endl << "Attention: the program only runs for the equilibration phase";
	}
	folder = "Data/" + state + "/" + distribution + "/";
	
	cout << endl << endl << "Starting position (in Bohr radius units): " << endl;
	cout << "x: " << x_0 << endl;
	cout << "y: " << y_0 << endl;
	cout << "z: " << z_0 << endl;
	cout << endl << "Delta: " << delta << endl;
	cout << "Atomic state: " << state << endl;
	cout << "Metropolis' distribution: " << distribution << endl; 
	

// *************** EQUILIBRATION ****************** 

	int N=100;	//per equilibrare faccio fare N step a vuoto
	int accept=0;
	
	for(int i=0; i<N; i++)	pos = Step(pos, &accept);
	
	cout << endl << "*** EQUILIBRATION PHASE ***" << endl;
	cout << endl << "The step has been accepted " << accept << " times over " << N << endl;
	cout << "The final position is: (" << pos.x << ", " << pos.y << ", " << pos.z << ")" << endl << endl;
	
	if(!start){
		ofstream equilibration(folder + "equilibration.dat", ios::app);	//ad ogni esecuzione del codice aggiungo r iniziale, delta e accettazione
		equilibration << R(x_0,y_0,z_0) << " " << delta << " " << accept/(double)N << endl;
		equilibration.close();
	}	
	//stampo solo se start == false 
	
// *************** <r> + position ******************

if(start){
	ofstream out_rad, out_pos;
	out_rad.open(folder + "radius.dat");	//raggio medio con errore all'avanzare dei blocchi
	out_pos.open(folder + "position.dat");	//posizione (x,y,z)
	
	int N_blocchi = 100;
	int M = 1000000;
	int L = M/N_blocchi;																				//numero di passi per ogni blocco
	
	double average=0, average2=0, error=0;
	for(int n=0; n<N_blocchi; n++){
		int accept=0;																							//numero di valori accettati per ogni blocco (per controllo)
		double block=0, block2=0;
		for(int istep=0; istep < L; istep++){
			pos = Step(pos, &accept);
			if(istep%10==0)	out_pos << pos.x << " " << pos.y << " " << pos.z << endl;	//stampo posizione ogni 10 passi
			block += R(pos.x, pos.y, pos.z);
		}
		block /= (double)L;
		block2 = block*block;
		
		average = (double)(((average)*n + block)/(n+1.));
  	average2 = (double)(((average2)*n + block2)/(n+1.));				//medie all'avanzare dei blocchi
  	if( n!=0)	error = sqrt( (average2 - average*average)/(double)n);
		
		out_rad << n+1 << " " << average << " " << error << endl;
	}

	out_rad.close();
	out_pos.close();
}

return 0;
}

// ******************** FUNZIONI *********************

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

// ******************* //

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
	cerr << endl << "ERROR: wrong state or distribution" << endl;		//in caso di errore la posizione viene riportata al centro
	position err_pos = {0, 0, 0};
	return err_pos; 
}

position Step_1s_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd.Uniform(-delta, delta);
	new_pos.y += rnd.Uniform(-delta, delta);
	new_pos.z += rnd.Uniform(-delta, delta);
	
	double A = min(1., state_1s(new_pos.x,new_pos.y,new_pos.z)/state_1s(pos.x,pos.y,pos.z));
	double m = rnd.Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_1s_gaus(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd.Gauss(0, delta/2.);
	new_pos.y += rnd.Gauss(0, delta/2.);
	new_pos.z += rnd.Gauss(0, delta/2.);
	
	double A = min(1., state_1s(new_pos.x,new_pos.y,new_pos.z)/state_1s(pos.x,pos.y,pos.z));
	double m = rnd.Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_2p_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd.Uniform(-delta, delta);
	new_pos.y += rnd.Uniform(-delta, delta);
	new_pos.z += rnd.Uniform(-delta, delta);
	
	double A = min(1., state_2p(new_pos.x,new_pos.y,new_pos.z)/state_2p(pos.x,pos.y,pos.z));
	double m = rnd.Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_2p_gaus(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd.Gauss(0, delta/2.);
	new_pos.y += rnd.Gauss(0, delta/2.);
	new_pos.z += rnd.Gauss(0, delta/2.);
	
	double A = min(1., state_2p(new_pos.x,new_pos.y,new_pos.z)/state_2p(pos.x,pos.y,pos.z));
	double m = rnd.Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_3d_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd.Uniform(-delta, delta);
	new_pos.y += rnd.Uniform(-delta, delta);
	new_pos.z += rnd.Uniform(-delta, delta);
	
	double A = min(1., state_3d(new_pos.x,new_pos.y,new_pos.z)/state_3d(pos.x,pos.y,pos.z));
	double m = rnd.Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

position Step_4f_unif(position pos, int*accept){
	position new_pos = pos;
	new_pos.x += rnd.Uniform(-delta, delta);
	new_pos.y += rnd.Uniform(-delta, delta);
	new_pos.z += rnd.Uniform(-delta, delta);
	
	double A = min(1., state_4f(new_pos.x,new_pos.y,new_pos.z)/state_4f(pos.x,pos.y,pos.z));
	double m = rnd.Rannyu();
	
	if(m<A){
		pos = new_pos;
		(*accept)++;
	}	
	return pos;
}

// ******************** //

double R(double x, double y, double z){
	return sqrt( x*x + y*y + z*z );
}

double cos_z(double x, double y, double z){
	return z/R(x,y,z);
}
