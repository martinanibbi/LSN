#include<vector>
#include "random.h"

using namespace std;

#ifndef _Travelling_Salesman_Simulated_Annealing_
#define _Travelling_Salesman_Simulated_Annealing_

//NB: P_Mutation inizializzato sempre a 0 e Mutation funziona in modo diverso rispetto alla lezione_09 
//----> inidividuo subisce PER FORZA una mutazione, scelta a caso con la stessa probabilità tra local_shift, permutation, swap

	//coordinate di ciascuna città estratte lungo una circ o uniformemente in un quadrato
struct coordinates{	
	double x;
	double y;
};

	//probabilità di tutte le mutazioni e di crossover 
struct probability{
	double crossover;
	double swap;
	double global_shift;
	double local_shift;
	double permutation;
	double inversion;
};

// ******************** //

class Individual{
	private:
		const int n_cities = 30;	
			//numero delle città -> 30
		vector<int> journey;	
			//ordine in cui le città vengono visitate		//ATTENZIONE: da 0 a 29 (indici che uso --> city[i]) (!!!!!!)
		vector<coordinates> city;	
			//coordinate delle città (devono essere le stesse per tutti gli individui!)
		probability P_Mutation;
			//struttura che contiene tutte le probabilità di mutazione e crossover (uguale per tutti gli individui!!)
		
		bool control;
		double path_length;
		Random* rnd;
			//Nota: tutti gli individui devono puntare allo stesso generatore random!
		
		double Distance(int step);	
			//"distanza" di due città contigue
		void Total_distance();			
			//misura L2 (quadrato della distanza) e salva in path_length
		
		void Forward_Global_Shift(int);
		void Backward_Global_Shift(int);
		void Forward_Local_Shift(int,int);
		void Backward_Local_Shift(int,int);		
			//utilizzati in Random ... Shift()
		
	public:
		Individual();
		Individual(const Individual&);
		Individual(Random*, vector<coordinates>, probability);
		Individual& operator=(const Individual&);
			//Nota: journey, city, rnd, P_Mutation
		void Check();								
			//se vincoli rispettati: control = true, altrimenti: control = false
		
		void Print(string);
			//stampo configurazione (devo fornire percorso del file!)
		
		void Set_Journey(vector<int> J){ journey = J; };
		void Set_City(vector<coordinates> c){ city = c; };
		vector<int> Get_Journey(){ return journey; };
		Random* Get_Rnd(){ return rnd; };
		vector<coordinates> Get_City(){ return city; }; 
		
		double Path_Length();				
			//eseguo Total_distance e restituisce il valore di path_length
		
		void Mutation();						
			//eseguo una mutazine tra swap, permutation e local shift
		void Swap(int,int);
			//scambia due elementi di journey
		void Random_Global_Shift();				
			//Shift non modifica il percorso ma serve per ergodicità!			
		void Random_Local_Shift();	
			//Local shift fa scalare m elementi con m tra [3,29]
		void Permutation();
			//Permutation scambia gruppi di m elementi (NB: m<n_cities/2)
		void Inversion();
			//Inversion inverte m elementi
			
		void Crossover(Individual parent1, Individual parent2);	
			//sostituisce journey locale t.c. Individual figlio con metà corredo genetico del genitore1 INVARIATO e completato con il DNA del  
			//genitore2 (rispettando i vincoli...)
};

// ******************** //

class Sim_Ann{
	private:
		Random* rnd;
		
		const int n_cities = 30; 
		
		string city_conf, folder;
			//può essere 'circle' oppure 'square'
		probability P_Mutation;
			//struttura che contiene tutte le probabilità di mutazione (NB: crossover e global shift ==0) --> lette in Input.dat
		
		double temp;
			//campionamento dipende da temperatura, inizializzata e modificata da main
		
		double acceptance[3], n_steps[3];
			//NB: va campionato solo con la stessa temperatura!!!!		//Non esiste delta, è una mossa "secca"
			//0-> SWAP		1->LOCAL SHIFT		2->PERMUTATION
		vector<coordinates> city;				
			//coordinate delle città
		
		Individual System;	
			//sistema fisico simulato con energia pari alla lunghezza del percorso
		Individual Best;
			//ad ogni esecuzione di Step() controllo se o ottenuto na configurazione migliore di Best
	
	public:
		Sim_Ann(Random*,double);
		Sim_Ann(Random*,double,string);	
			//ricavo configurazione città da file ----> RICAVATE DA ES_01!!!!
		Sim_Ann(Sim_Ann&);
			//copy constructor
		
		void Set_Temp(double t);
			//NB: ogni volta che cambio la temperatura devo resettare acceptance e n_steps
		void Set_City(vector<coordinates> c){ System.Set_City(c); };
		string Get_City_Conf(){ return city_conf; };
			//per sapere se mi trovo in 'square' o 'circle'
		vector<coordinates> Get_City(){ return System.Get_City(); };
		
		void Input();
			//chiamata dal costruttore, legge da input.dat i dati preliminari
		void City_constructor();
			//chiamata dal costruttore, genera configurazione di città in un quadrato o su una circonferenza (decide in Input), non stampo
		void City_constructor(string);
			//ricavo configurazione città da file (DA ESERCIZIO 1), non stampo
			
		void Step();
			//eseguo uno step con Metropolis (ad una temperatura temp fissata)
		int Mutation(Individual&);
			//di fatto è un overload di Mutation di Individual ---> eseguo SEMPRE una mutazione, solo scelta a caso tra swap, local shift e 
			//permutation (con la stessa percentuale...)
			//restituisco 0 per swap, 1 per local_shift, 2 per permutation ---> salvo valori di accettazione in 3 variabili diverse
			
		double Path_Length(){	return System.Path_Length(); };	
			//restituisco lunghezza dell'n-esimo percorso
		double Best_Path_Length(){ return Best.Path_Length(); };
			//restituisco lunghezza del miglior percorso visitato
		void Print_City();
		void Print_Path();
		void Print_Path(int size, int rank);
			//stampo ricorsivamente percorso alla temperatura T ---> in file diversi per ogni rank!
		void Print_Conf();
		void Print_Conf(int size, int rank);
			//stampo configrazione del percorso corrente --> in file diversi per ogni rank!
};

#endif






