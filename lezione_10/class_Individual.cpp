#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
#include "Travelling_Salesman.h"

using namespace std;

// *********************************************** Inizializzazione, controllo, stampa ***********************************************

Individual :: Individual(){
	rnd = NULL;
	P_Mutation = {0,0,0,0,0,0};
}	

	//inizializza vettore con configurazione casuale	//NB: alla fine devo controllare con Check()
	//Questo metodo è più ergodico rispetto allo scambio di tot elementi partendo da una configurazione ordinata
	
Individual :: Individual(Random* r, vector<coordinates> c, probability p){
	rnd = r;	//rnd presente in population... tutti gli individui puntano allo stesso generatore!!
	city = c;	//copio in ogni elemento la configurazione delle città (proviene da population ed è uguale per ogni individuo)
	P_Mutation = p;	//idem x probabilità di mutazioni
	
	vector<int> available_cities;		//Carico una configurazione delle città ordinata da 0 a n_cities-1
	for(int i=0; i<n_cities; ++i){
		available_cities.push_back(i);
	}	
	for(int i=0; i<n_cities; ++i){				//carico il vettore journey
		int n_max = n_cities - i;
		double x = rnd->Rannyu(0., (double)n_max);					//estraggo numero tra 0 e (n_cities-i) ---> al primo passaggio (0,30), (0,29)...
		int index = (int)x;																//indice di available_cities è parte intera del numero estratto	(va da 0 a nmax-1)	
		journey.push_back( available_cities[index] );			//assegno la città
		available_cities.erase(available_cities.begin() + index);				
				//ad ogni passaggio cancello l'elemento usato e riduco la dimensione di available_cities	
	}		
	Check();	//controllo se non si ripetono elementi (control = true se tutto ok)
	Total_distance();		//salvo lunghezza iniziale del percorso in path_length
}

//*******************

	//serve per scambio in Population...
Individual& Individual :: operator = (const Individual& I){
	journey = I.journey;
	rnd = I.rnd;
	city = I.city;
	P_Mutation = I.P_Mutation;
	//Total_distance();	//calcolo e salvo nuova distanza
	return *this;
}

Individual :: Individual(const Individual& I){
	journey = I.journey;
	rnd = I.rnd;
	city = I.city;
	P_Mutation = I.P_Mutation;
}


//*******************

void Individual :: Check(){	
	for(int i=0; i<n_cities-1; ++i){
		for(int j=i+1; j<n_cities; ++j){
			if(journey[i] == journey[j]){
				cerr << endl << "Error: journey vector out of bounds" << endl;
				control = false;
				return;
			}
		}
	}
control = true;	
return;
}

//******************

void Individual :: Print(string filename){
	ofstream OutJ(filename);
	for(int i=0; i<n_cities; ++i)	OutJ << journey.at(i) << endl;
	OutJ.close();
return;
}

// *********************************************** Lunghezza del percorso ***********************************************

double Individual :: Path_Length(){	
	Total_distance();
	return path_length;
}

void Individual :: Total_distance(){	//NOTA: è il quadrato della lnghezza!
	double dist = 0.;
	for(int step=1; step<=n_cities; step++){
		dist += Distance(step);
	}
path_length = dist;	
}

	//all'i-esimo step: | ( x_i-x_(i-1) )^2 + ( y_i-y_(i-1) )^2 |
	//Attenzione: all'ultimo step devo tornare all'inizio!
double Individual :: Distance(int step){
	double dist = 0.;
	if(step == n_cities){
		dist += ( (city.at(journey[n_cities-1])).x - (city.at(journey[0])).x )*( (city.at(journey[n_cities-1])).x - (city.at(journey[0])).x );
		dist += ( (city.at(journey[n_cities-1])).y - (city.at(journey[0])).y )*( (city.at(journey[n_cities-1])).y - (city.at(journey[0])).y );
	}
	else{
		dist += ( (city.at(journey[step])).x - (city.at(journey[step-1])).x )*( (city.at(journey[step])).x - (city.at(journey[step-1])).x );
		dist += ( (city.at(journey[step])).y - (city.at(journey[step-1])).y )*( (city.at(journey[step])).y - (city.at(journey[step-1])).y );
	}
return dist;
}

// *********************************************** Mutazioni di journey ***********************************************

void Individual :: Swap(int i, int j){		//scambio elementi i - j di journey
	i = i%n_cities;
	j = j%n_cities;	//applico le periodic boundary conditions
	//cout << i << "	" << j << endl;
	if(i == j)	return;

	int appo;
	appo = journey.at(i);
	journey.at(i) = journey.at(j);
	journey.at(j) = appo;
	return;
}

// ****************** //

void Individual :: Random_Global_Shift(){	//posso scalare da 1 a 15, in avanti o indietro
	int r = rnd->Rannyu(1., n_cities/2.+1.);
	double s = rnd->Rannyu(-1.,1.);
	if(s > 0)	Forward_Global_Shift(r);	
	else			Backward_Global_Shift(r);
	//cout << endl << "Mutation: " << r << "-shift" << endl;
	return;
}

void Individual :: Forward_Global_Shift(int N){	//scalo di N posti in avanti		
	for(int n=0; n<N; n++){	
		for(int i=0; i<n_cities-1; i++)		Swap(i,i+1);
	}
	return;
}

void Individual :: Backward_Global_Shift(int N){	//scalo di N posti all'indietro	
	for(int n=0; n<N; n++){	
		for(int i=n_cities-1; i>1; i--)		Swap(i,i-1);
	}
	return;
}

// **************** //

void Individual :: Random_Local_Shift(){		//scelgo porzione da far scalare (3,29), numero di passi (1, n/2 + 1), avanti o indietro
	int m = rnd->Rannyu(3,30);
	int r = rnd->Rannyu(1., (double)m/2.+1.);
	double s = rnd->Rannyu(-1.,1.);
	if(s > 0)	Forward_Local_Shift(m,r);	
	else			Backward_Local_Shift(m,r);
	return;
}

void Individual :: Forward_Local_Shift(int m, int N){		//faccio scalare una porzione di m elementi di N posti in avanti		
	int begin = rnd->Rannyu(0,30);	//dove incomincio a scambiare
	int end = begin + m;						//dove finisco
	for(int n=0; n<N; n++){	
		for(int i=begin; i<end; i++)		Swap(i,i+1);
	}
	return;
}

void Individual :: Backward_Local_Shift(int m, int N){	//faccio scalare una porzione di m elementi di N posti all'indietro
	int begin = rnd->Rannyu(0,30);	//dove incomincio a scambiare
	int end = begin + m;						//dove finisco
	for(int n=0; n<N; n++){	
		for(int i=end; i>begin; i--)		Swap(i,i-1);
	}
	return;
}

// **************** //

void Individual :: Permutation(){	//permuta m elementi contigui		//NB: non devono esserci sovrapposizioni --> m<dim_population/2
	int m = rnd->Rannyu(2, (double)n_cities/2.);	// m -> numero elementi da permutare con altri m elementi contigui
	int begin = rnd->Rannyu(0,30);	//dove incomincio a scambiare
	int end = begin + m;
	for(int i=begin; i<end; ++i){
		Swap(i,i+m);
	}
	return;
}

// **************** //

void Individual :: Inversion(){
	int m = rnd->Rannyu(2, (double)n_cities/2.);	// m -> numero elementi da invertire con altri m elementi contigui
	int begin = rnd->Rannyu(0,30);	//dove incomincio a invertire (indice da 0 a 29)
	int end = begin + m;
	for(int i=begin; i<end/2; ++i){
		Swap(i,end-i);
	}
	return;
}

// **************** //

void Individual :: Mutation(){	//scelgo quale/i mutazione/i fare avvenire
	
	double r_global_shift = rnd->Rannyu();			//ogni mutazione avviene con tot probabilità
	if( r_global_shift < P_Mutation.global_shift ){								
		Random_Global_Shift();
	}	
	
	double r_swap = rnd->Rannyu();		//Swap è la mutazione che contribuisce di più all'abbassamento di path_length
	if( r_swap < P_Mutation.swap){
		int sx = rnd->Rannyu(0.,30.);
		int dx = rnd->Rannyu(0.,30.);
		Swap(sx, dx);
	}
	
	double r_local_shift = rnd->Rannyu();
	if ( r_local_shift < P_Mutation.local_shift ){
		Random_Local_Shift();
	}
	
	double r_permutation = rnd->Rannyu();
	if ( r_permutation < P_Mutation.permutation ){
		Permutation();
	}
	
	double r_inversion = rnd->Rannyu();
	if ( r_inversion < P_Mutation.inversion){
		Inversion();
	}
	Check();
	return;
}

// *********************************************** Crossover ***********************************************

void Individual :: Crossover(Individual parent1, Individual parent2){	//NB: modifico journey locale!!
	double r = rnd->Rannyu();
	if(r<P_Mutation.crossover)	return;
		//crossover non avviene sempre!!!!
	
	vector<int> son = parent1.Get_Journey();	//inizialmente figlio = genitore1
	int m = n_cities/2 + 1;		//posizione in cui taglio
	
	for(int i=n_cities-1; i>=m; --i)	son.erase(son.begin()+i);	//cancello ultimi elementi fino al taglio
	
	vector<int> available_genes = parent2.Get_Journey();	//copio geni parent2, elimino geni già presenti in son
	for(int j=0; j<son.size(); ++j){	//ciclo su elementi di son
		for(int i=0; i<available_genes.size(); ++i){	//ciclo su elementi di available_genes
			if(available_genes.at(i) == son.at(j)){
				available_genes.erase(available_genes.begin()+i);
				break;
			}
		}
	}
	for(int i=0; i<available_genes.size(); ++i)	son.push_back(available_genes.at(i));
	journey = son;	
	return;
}	


