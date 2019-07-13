#include<vector>
#include "random.h"

using namespace std;

#ifndef _Travelling_Salesman_Genetic_
#define _Travelling_Salesman_Genetic_

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
		vector<int> Get_Journey(){ return journey; };
		Random* Get_Rnd(){ return rnd; };
		vector<coordinates> Get_City(){ return city; }; 
		
		double Path_Length();				
			//eseguo Total_distance e restituisce il valore di path_length
		
		void Mutation();						
			//chiamato in Population :: Breeding(), scambia elementi del vettore journey
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

class Population{
	private:
		Random* rnd;
		const int n_cities = 30;					
			//numero delle città -> 30
		const int dim_population = 900;		
			//popolazione -> n_cities^2
		int p, n_steps;									
			//potenza utilizzata nella scelta del genitore e numero di iterazioni in New_Generation()
		string city_conf, folder;
			//può essere 'circle' oppure 'square'
		probability P_Mutation;
			//struttura che contiene tutte le probabilità di mutazione e crossover (uguale per tutti gli individui!!) --> lette in Input.dat
		
		vector<coordinates> city;				
			//coordinate delle città (devono essere le stesse per tutti gli individui!)
		vector<Individual> single;			
			//vettore di individui
		Individual Best_single;
			//salvo il miglior elemento
		
		void Swap(int, int);
			//scambio due individui nel vettore single (NB: non più ordinato...)
		void Quick_Sort(int, int);		
			//usato in Quick_Sort() --> pubblico
		void Scambia_Min(int);				
			//usato in Selection_Sort() --> pubblico
		
	public:
		Population(Random*);
		void Input();
			//chiamata dal costruttore, legge da input.dat i dati preliminari
		void City_constructor();
			//chiamata dal costruttore, genera configurazione di città in un quadrato o su una circonferenza (decide in Input), stampo 
		
		double Path_Length(int n){	return single[n].Path_Length(); };	
			//restituisco lunghezza dell'n-esimo percorso
		
		void Quick_Sort();				
			//ordina single in base alla distanza
		void Selection_Sort();
		
		void Evolution();
			//eseguo New_generation per n_steps volte
		void New_generation();		
			//rimpiazzo tutto single
		void Breeding(vector<Individual>&);							
			//ordino, scelgo genitori, elimino elementi peggiori, scambio/muto DNA, aggiungo figli al vettore new_gen
		Individual Parent_choice();		
			//scelgo genitori che fittano meglio ---> NB: devo ordinare il vettore FUORI da questo metodo! 
					
		void Print_Average_Path(int);
		void Print_Best_Path(int);
			//stampano (ricorsivamente) media della migliore metà e miglior percorso 
			//Nota: ricevono in ingresso il "passo" dell'evoluzione											
};

#endif
