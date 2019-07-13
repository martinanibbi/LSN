#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
#include "Travelling_Salesman.h"

using namespace std;

// *********************************************** Inizializzazione Population ***********************************************

Population :: Population(Random* r){

	rnd = r;
	Input();	//leggo configurazione delle città e probabilità di mutazioni/crossover
	City_constructor();	//costruisco configurazione delle città in un/a quadrato/circonferenza	
	
	for(int i=0; i<dim_population; ++i){
		Individual ind(rnd,city,P_Mutation);	//inizializzo e riempio vector di 900 individui
		single.push_back(ind);
	}
	Best_single = single[0];	//inizializzato a caso, lo aggiorno ogni volta che eseguo Quick_Sort()
	Quick_Sort();
}

// ******************* //

void Population :: Input(){
	ifstream ReadInput("Data/input.dat");
		ReadInput >> city_conf;
			folder = "Data/" + city_conf + "/";
		ReadInput >> P_Mutation.crossover;
		ReadInput >> P_Mutation.swap;
		ReadInput >> P_Mutation.global_shift;
		ReadInput >> P_Mutation.local_shift;
		ReadInput >> P_Mutation.permutation;
		ReadInput >> P_Mutation.inversion;
		ReadInput >> p;
		ReadInput >> n_steps;
	ReadInput.close();
	return;
}

// ******************* //

void Population :: City_constructor(){
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
	ofstream OutCity(folder + "city_conf.0");
	for(int i=0; i<n_cities; ++i)	OutCity << city.at(i).x << " " << city.at(i).y << endl;
	
	return;
}

// *********************************************** Ordinamento Population ***********************************************

void Population :: Selection_Sort(){
 	for(int i=0; i<dim_population-1; i++){
		Scambia_Min(i);
	}	
	return;
}

void Population :: Scambia_Min(int i){

double min = single[i].Path_Length();
int iMin=i;

	for(int j=i+1; j<dim_population; j++){
		if(single[j].Path_Length() < min){
			min = single[j].Path_Length();
			iMin=j;
		}
	}

Swap(i,iMin);
}

// ******************* //

void Population :: Quick_Sort(){
	Quick_Sort(0, dim_population-1);	
	if(single[0].Path_Length() < Best_single.Path_Length() )	Best_single = single[0];
	return;
}

void Population :: Quick_Sort(int primo, int ultimo){

	if(primo>ultimo || ultimo>=dim_population) return;	// controllo validita' degli estremi
	if( ultimo - primo <= 1 ){
		if( single[primo].Path_Length() > single[ultimo].Path_Length() ) Swap(primo,ultimo);
		return;
	}

	int index = int(primo+ultimo)/2;	//indice dell'elemento in mezzo al vector single
	double pivot = single[ index ].Path_Length();	//path_length dell'elemento in mezzo al vettore
	int basso = primo;
	int alto = ultimo;
	
	while(basso<alto){	//tutti gli elementi minori del pivot vanno a sx e viceversa
		
		while( (single[basso].Path_Length()) < pivot)	basso++;	
				// aumenta indice basso fino a quando lunghezze sono piu' piccole rispetto a pivot 
		while( (single[alto].Path_Length()) > pivot)  alto--;	
				// diminuisce indice alto fino a quando lunghezze sono piu' grandi rispetto a pivot
		
		// ora ho le prime componenti "fuori posto" con indici alto e basso
		if(basso<alto){
			Swap(basso,alto);	
			basso++;
		}
	}
	// ripeto il tutto, divido due intervalli primo - basso-1, basso - alto
	Quick_Sort(primo,basso-1) ;
	Quick_Sort(basso,ultimo);
}

// **************** //
	
void Population :: Swap(int a, int b){	//scambio individui del vettore single
	Individual appo(single.at(a));
	single.at(a) = single.at(b);
	single.at(b) = appo;
	return;
}

// *********************************************** Evoluzione ***********************************************

void Population :: Evolution(){
	for(int i=1; i<=n_steps; ++i){
		New_generation();	//vettore già ordinato
		if(i%10 ==0){
			Print_Average_Path(i);
			Print_Best_Path(i);	//metodi Print NON eseguono Quick_Sort(), già ordinato in New_generation() 
		}
	}
	string filename = folder + "best_path_conf.dat";
	single.at(0).Print(filename);
		//stampo miglior configurazione 
	for(int j=0; j<dim_population; ++j) cout << (single.at(j)).Path_Length() << endl;
		//stampo a video per l'ultimo passaggio...	
	return;
}

// **************** //

void Population :: New_generation(){
	vector<Individual> new_gen;
	
	for(int i=0; i<dim_population/2.; i++){
		Breeding(new_gen);	
			//ad ogni passaggio aggiungo 2 elementi a new_gen, fino a sostituire single
	}
	single = new_gen;
	Quick_Sort();
	single.erase(single.begin()+dim_population-1);	
	single.insert(single.begin(), Best_single);	//sostituisco il migliore al peggiore, così single rimane ordinato
	return;
}

// **************** //

void Population :: Breeding(vector<Individual>& new_gen){
	//Quick_Sort();	//ordino single fuori da questo metodo
	
	Individual father = Parent_choice();
	Individual mother = Parent_choice();
		
	Individual son1(father);
		son1.Crossover(father,mother);
	Individual son2(mother); 
		son2.Crossover(mother,father);	
		
			//non importa come sono inizializzati son1 e son2, journey viene assegnato definitivamente in Crossover
			//NB: genitore 1 e genitore2 hanno due ruoli leggermente diversi... 
	
	//NB: sostituisco gli ultimi elementi (potrebbero anche essere i genitori...) del vettore single con i figli 
	son1.Mutation();
	son2.Mutation();		//non sempre avviene mutazione...
	
	new_gen.push_back(son1);
	new_gen.push_back(son2);
		//NB: riempio vettore new_gen
	return;
}

// **************** //

Individual Population :: Parent_choice(){				//NB: devo ordinare il vettore FUORI da questo metodo!
	double r = rnd->Rannyu();
	int j = int(single.size()*(pow(r,p)));				//r compreso tra 0 e 1, più probabile che il numero sia piccolo (-->path_length minore)
																								//NB: dimensione di single potrebbe cambiare!!!
	return single.at(j); 	
}

// *********************************************** Risultati ***********************************************

void Population :: Print_Best_Path(int step){
	//Quick_Sort();
	ofstream OutPath(folder + "best_path_length.dat", ios::app);
	OutPath << step << " " << single.at(0).Path_Length() << endl;
	OutPath.close();	
	return;
}

// **************** //

void Population :: Print_Average_Path(int step){
	//Quick_Sort();
	double sum = 0.;
	for(int i=0; i<dim_population; ++i){
		sum += single.at(i).Path_Length();
	}
	sum /= (double)(dim_population)/2.;
	ofstream OutAve(folder + "average_path_length.dat", ios::app);
	OutAve << step << " " << sum << endl;
	OutAve.close();
	return;
}

