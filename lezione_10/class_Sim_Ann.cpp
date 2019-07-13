#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
#include "Travelling_Salesman.h"

using namespace std;

// *********************************************** Inizializzazione System ***********************************************

Sim_Ann :: Sim_Ann(Random* r,double t){
	temp = t;
	rnd = r;
	Input();	//leggo configurazione delle città e probabilità di mutazioni/crossover
	City_constructor();	//costruisco configurazione delle città in un/a quadrato/circonferenza	
	
	for(int i=0; i<3; ++i){
		acceptance[i] = 0;
		n_steps[i] = 0;
	}	
	Individual ind(rnd,city,P_Mutation);	//inizializzo System (a caso)
	System = ind;
	Best = ind;
}

// ******************* //

Sim_Ann :: Sim_Ann(Random* r, double t, string filename){
	temp = t;
	rnd = r;
	Input();
	City_constructor(filename);	
	
	for(int i=0; i<3; ++i){
		acceptance[i] = 0;
		n_steps[i] = 0;
	}	
	Individual ind(rnd,city,P_Mutation);	//inizializzo System (a caso)
	System = ind;
	Best = ind;
}

// ******************* //

Sim_Ann :: Sim_Ann(Sim_Ann& copy){
	temp = copy.temp;
	rnd = copy.rnd;
	P_Mutation = copy.P_Mutation;
	city_conf = copy.city_conf;
	folder = copy.folder;
	city = copy.city;
	System = copy.System;
	Best = copy.Best;
	for(int i=0; i<3; ++i){
		acceptance[i] = 0;
		n_steps[i] = 0;
	}	
}

// ******************* //

void Sim_Ann :: Input(){
	ifstream ReadInput("Data/input.dat");
		ReadInput >> city_conf;
			folder = "Data/" + city_conf + "/";
		ReadInput >> P_Mutation.crossover;
		ReadInput >> P_Mutation.swap;
		ReadInput >> P_Mutation.global_shift;
		ReadInput >> P_Mutation.local_shift;
		ReadInput >> P_Mutation.permutation;
		ReadInput >> P_Mutation.inversion;
	ReadInput.close();
	return;
}

// ******************* //

void Sim_Ann :: City_constructor(){
	if(city_conf == "square"){
		for(int i=0; i<n_cities; ++i){
			double X = rnd->Rannyu(-1,1);
			double Y = rnd->Rannyu(-1,1);
			coordinates position = {X,Y};
			city.push_back(position);
		}
	}
	if(city_conf == "circle"){
		for(int i=0; i<n_cities; ++i){
			double theta = rnd->Rannyu(0,2*M_PI);	//estraggo theta uniformemente tra 0 e 2pi
			coordinates position = {cos(theta), sin(theta)};
			city.push_back(position);
		}
	}	
	
	return;
}

// ******************* //

void Sim_Ann :: City_constructor(string filename){
	ifstream ReadCity(filename);
	if(ReadCity.fail())	cerr << endl << "Error: unable to read city configuration" << endl;
	city.clear();
	double X,Y;
	for(int i=0; i<n_cities; ++i){
		ReadCity >> X;
		ReadCity >> Y;
		coordinates position = {X,Y};
		city.push_back(position);
	}
	
	ReadCity.close();
	
	return;
}

// ******************* //

void Sim_Ann :: Set_Temp(double t){
	temp = t;
	for(int i=0; i<3; i++){
		acceptance[i] = 0;
		n_steps[i] = 0;
	}
return;
}

// *********************************************** Metropolis ***********************************************

void Sim_Ann :: Step(){
	Individual appo(System);
	int index = Mutation(appo);
	double Delta_E = appo.Path_Length() - Path_Length();	//energia = lunghezza del percorso nuovo - lunghezza del percorso vecchio 
	double beta = 1./temp;
	double A = min( 1., exp(-beta*Delta_E) );
	double m=rnd->Rannyu();
	if(m<A){
		System = appo;
		acceptance[index]+=1.;
		if(System.Path_Length() < Best.Path_Length())	Best = System;
	}
	n_steps[index]+=1.;
	
	
	return;
}

// ******************* //

int Sim_Ann :: Mutation(Individual& appo){	//scelgo quale mutazione fare avvenire sull'elemento ausiliario appo
	
	int pippo = rnd->Rannyu(0,3);			//ogni mutazione avviene con probabilità di 1/3
	
	if(pippo == 0){
		int sx = rnd->Rannyu(0.,30.);
		int dx = rnd->Rannyu(0.,30.);
		appo.Swap(sx, dx);
	}
	if (pippo == 1){
		appo.Random_Local_Shift();
	}
	if (pippo == 2){
		appo.Permutation();
	}
	appo.Check();
	return pippo;
}

// *********************************************** Risultati ***********************************************
	
	//NB: stampo il migliore!!!!!
void Sim_Ann :: Print_Path(){
	//Quick_Sort();
	ofstream OutPath(folder + "path_length.dat", ios::app);
	OutPath << temp << " " << Best.Path_Length() << endl;
	OutPath.close();
	vector<double> acc_per;
		for(int i=0; i<3; i++)	acc_per.push_back(100.*acceptance[i]/(double)n_steps[i]);
		
	cout << endl << "Temperatura:	" << temp << endl;
	cout << "Accettazione swap:		" << acc_per.at(0) << endl;
	cout << "Accettazine local shift:	" << acc_per.at(1) << endl;
	cout << "Accettazione permutazione:	" << acc_per.at(2) << endl;
	cout << "Lunghezza del percorso:	" << Path_Length() << endl;
	return;
}

//**********************

void Sim_Ann :: Print_City(){
	ofstream OutCity(folder + "city_conf.0");
	for(int i=0; i<n_cities; ++i)	OutCity << city.at(i).x << " " << city.at(i).y << endl;
	OutCity.close();
	return;
}

//**********************

void Sim_Ann :: Print_Path(int size, int rank){
	//Quick_Sort();
	for(int i=0; i<size; i++){
		if(i==rank){
			ofstream OutPath(folder + "path_length_" + to_string(rank) + ".dat", ios::app);
			OutPath << temp << " " << Best.Path_Length() << endl;
			OutPath.close();
		}
	}
	/*vector<double> acc_per;
		for(int i=0; i<3; i++)	acc_per.push_back(100.*acceptance[i]/(double)n_steps[i]);
	
	
	cout << endl << "Temperatura:	" << temp << endl;
	cout << "Accettazione swap:		" << acc_per.at(0) << endl;
	cout << "Accettazine local shift:	" << acc_per.at(1) << endl;
	cout << "Accettazione permutazione:	" << acc_per.at(2) << endl;
	cout << "Lunghezza del percorso:	" << Path_Length() << endl;
	*/
	return;
}

//**********************

void Sim_Ann :: Print_Conf(){
	string filename = folder + "best_path_conf.dat";
	Best.Print(filename);
		//stampo miglior configurazione 
}

//**********************

void Sim_Ann :: Print_Conf(int size, int rank){
	for(int i=0; i<size; i++){
		if(i==rank){
			string filename = folder + "best_path_conf_" + to_string(rank) + ".dat";
			Best.Print(filename);
				//stampo miglior configurazione 
		}
	}
}
