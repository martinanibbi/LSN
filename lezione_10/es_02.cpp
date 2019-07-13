#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include "random.h"
#include "Travelling_Salesman.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]){
	MPI::Init(argc,argv);
		
		int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();
		
		// ************************ Inizializzazione Random **********************
		
		Random rnd;
		int seed[4];
	  int p1, p2;
	  for(int i=0; i<size; i++){	//ogni rank legge 2 primes diversi (senza rischio di apertura contemporanea di file)
		  if(i==rank){
 				ifstream Primes("Primes");
			  if (Primes.is_open()){
			    for(int j=0; j<rank; j++)	Primes >> p1 >> p2 ;
			  } else cerr << "PROBLEM: Unable to open Primes" << endl;
			  Primes.close();
			}
		}

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

		// *********************** Inizializzazione Simulated Annealing ************************
		
		double* recv = new double[size];
		for(int i=0; i<size; i++){
			recv[i] = 0;
		}
		
		Sim_Ann SA(&rnd, 10);
		string conf = SA.Get_City_Conf();
	
		//NB: di default il costruttore di Simm_Ann genera configurazione casuale di città... ma primes diversi per ogni core quindi città diverse!
		//rank 0 legge configurazione di città generata nell'esercizio 1 e Bcast agli altri rank
		//siccome la informazione completa è racchiusa in un vector di strutture, inizializzo una struttura di coordinate in MPI
	
		const int nitems=2;
    int          blocklengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION};
    MPI_Datatype mpi_coor_type;
    MPI_Aint     offsets[2];

    offsets[0] = offsetof(coordinates, x);
    offsets[1] = offsetof(coordinates, y);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_coor_type);
    MPI_Type_commit(&mpi_coor_type);
		
		vector<coordinates> city;
		city.resize(30);
		
    if (rank == 0) {

        ifstream ReadCity("Data/" + conf + "/city_conf.0");
        for(int i=0; i<30; i++){
        	double X, Y;
        	ReadCity >> X >> Y;
        	coordinates pos = {X,Y};
        	city.at(i) = pos;
        }
        ReadCity.close();
        //for(int i=0; i<30; i++) cout << (city.at(i)).x << "	" << (city.at(i)).y << endl;
    }
    MPI_Bcast(&city.front(), 30, mpi_coor_type, 0, MPI_COMM_WORLD);
    SA.Set_City(city);
    
		
		// ******************************** Simulated Annealing ********************************
		
		for(double i=10; i>1; i-=0.05){
			SA.Set_Temp(i);
			for(int j=0; j<1000; j++)		SA.Step();
			SA.Print_Path(size, rank);	//stampo miglior elemento dopo 10000 step (per ogni rank si usano file diversi)
		}
		for(double i=1; i>0; i-=0.005){
			SA.Set_Temp(i);
			for(int j=0; j<10000; j++)	SA.Step();
			SA.Print_Path(size, rank);
		}
		for(double i=0.0049; i>0; i-=0.0001){
			SA.Set_Temp(i);
			for(int j=0; j<10000; j++)	SA.Step();
			SA.Print_Path(size, rank);
		}
		SA.Print_Conf(size, rank);
		double best = SA.Best_Path_Length();	//ottengo miglior percorso		
		
		// *******************
		
		for(int i=0; i<size; i++){
			if(i==rank){
				ofstream out("Data/" + conf + "/Parallel_Simulated_Annealing.dat", ios::app);
				out << rank << " " << best << endl;	//stampo rank e risultato finale		
				cout << rank << " " << best << endl;		
				out.close();
			}
		}
		
		MPI_Gather(&best,1,MPI_DOUBLE_PRECISION,recv,1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);
		double comm_best = recv[1];	//migliore tra tutti i percorsi ottenuti con calcolo parallelo
		
		if(rank == 0){
			for(int i=0; i<size; i++)	if(comm_best > recv[i])	comm_best = recv[i];		
			ofstream out("Data/" + SA.Get_City_Conf() + "/Parallel_Simulated_Annealing.dat", ios::app);
			out << size << " " << comm_best << endl;	//aggiungo al file il miglior risultato nella posizione 0		
			cout << endl << size << " " << comm_best << endl;
			out.close();
		}	
		
		delete[]recv;

  MPI_Type_free(&mpi_coor_type);
	MPI::Finalize();
	
return 0;
}
