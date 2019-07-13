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

		// ********************

	double command;
	if(rank == 0) command = 1.;
	else					command = 0.;
	
	MPI_Bcast(&command, 1, MPI_DOUBLE_PRECISION, 0, MPI::COMM_WORLD);
	
	for(int i=0; i<size; i++){
		if(i==rank)	cout << endl << i << "	" << command << endl;
	}
	
	
	
		const int nitems=2;
    int          blocklengths[2] = {1,1};
    MPI_Datatype types[2] = {MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION};
    MPI_Datatype mpi_coor_type;
    MPI_Aint     offsets[2];

    offsets[0] = offsetof(coordinates, x);
    offsets[1] = offsetof(coordinates, y);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_coor_type);
    MPI_Type_commit(&mpi_coor_type);

		const int tag = 16;
    if (rank == 0) {
        vector<coordinates> send;
        ifstream ReadCity("Data/square/city_conf.0");
        for(int i=0; i<30; i++){
        	double X, Y;
        	ReadCity >> X >> Y;
        	coordinates pos = {X,Y};
        	send.push_back(pos);
        }
        ReadCity.close();
        for(int i=0; i<30; i++) cout << (send.at(i)).x << "	" << (send.at(i)).y << endl;

        const int dest = 1;
        MPI_Send(&send.at(0),   30, mpi_coor_type, dest, tag, MPI_COMM_WORLD);

        cout << endl << "Rank " << rank << endl;
    }
    if (rank == 1) {
        MPI_Status status;
        const int src=0;

        vector<coordinates> recv;
        recv.resize(30);

        MPI_Recv(&recv.at(0), 30, mpi_coor_type, src, tag, MPI_COMM_WORLD, &status);
        cout << endl << "Rank 1 received:" << endl;
        //cout << endl << recv.size() << endl;
        
        for(int i=0; i<30; i++){
        	//cout << "Tutto bene!" << endl;
        	cout << (recv.at(i)).x << "	" << (recv.at(i)).y << endl;
        	        	//cout << "Tutto bene!" << endl;
        }
                
    }

    MPI_Type_free(&mpi_coor_type);
		MPI::Finalize();


return 0;
}

/*

Sim_Ann SA(&rnd,10);
		string filename = "Data/" + SA.Get_City_Conf() + "/city_conf.0";
		
		int command = 0;	
		int*itag = new int[size-1];
		for(int i=0; i<size-1; i++)	itag[i] = i;
		
		for(int i=0; i<size; i++){
			if(i == rank){
				if(i != 0) MPI::COMM_WORLD.Recv(&command, 1, MPI::INTEGER, i-1, itag[i-1]);
				SA.City_constructor(filename);
				if(i != size-1) MPI::COMM_WORLD.Send(&command, 1, MPI::INTEGER, i+1, itag[i]);
			}
		}
		
*/


/*
Sim_Ann SA(&rnd,10);
							//cout << endl << "Fin qui tutto bene! " << endl;	
		vector<coordinates> city(30);
		double cityx[30];
		double cityy[30];
		if(rank == 0){
			SA.City_constructor("Data/" + SA.Get_City_Conf() + "/city_conf.0");	
			city = SA.Get_City();
			for(int i=0; i<30; i++){
				cityx[i] = (city.at(i)).x;
				cityy[i] = (city.at(i)).y;
			}
		}
		MPI_Bcast(cityx,30,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);	//Broadcast da 0 verso gli altri
		MPI_Bcast(cityy,30,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);	
		
		for(int i=0; i<30; i++){
			(city.at(i)).x = cityx[i];
			(city.at(i)).y = cityy[i];		
		}
		SA.Set_City(city);	
							//ricostruisco città da esercizio 01 (altrimenti tutte configrazioni diverse) nel rank 0,
							//spezzetto e spedisco la configurazione tramite broadcast agli altri
							//ricostruisco e assegno configurazione città
							//(NB: non posso far leggere tutti dallo stesso file!!!)
							
*/
