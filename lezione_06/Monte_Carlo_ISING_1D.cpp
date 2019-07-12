/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"
#include "random.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration and acceptance rate (for Metropolis)

	ofstream out[4];
	if(h == 0)		out[iu].open("Data/" + metro + "/energy/final.dat" ,ios::app);
  if(h == 0) 		out[ic].open("Data/" + metro + "/heat_capacity/final.dat",ios::app);
  if(h == 0.02)	out[im].open("Data/" + metro + "/magnetization/final.dat",ios::app);
  if(h == 0)		out[ix].open("Data/" + metro + "/susceptibility/final.dat",ios::app);		//...altrimenti non apro file
    
  for(int i=0; i<4; i++){	//stampo solo risultati finali in funzione della temperatura
  	out[i] << temp << " " << glob_av[i]/(double)nblk << " " << err[i] << endl;
    out[i].close();
	}

  return 0;
}

// ************************ Input ************************

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("Data/input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // Metropolis or Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro == "Metropolis") 		cout << "The program perform Metropolis moves" << endl;
  else if(metro == "Gibbs")			cout << "The program perform Gibbs moves" << endl;
  else													cerr << "Wrong Input" << endl;
  
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  
  string start, eq;
  ReadInput >> start;		//leggo se devo ripartire da vecchia configurazione o nuova
  ReadInput >> eq;			//leggo se devo equilibrare (se parto da nuova conf di default...)
  
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
 	if(start == "old"){
 		Config_Old();
 		if(eq == "y")	Equilibration();
 	}
 	else{
 		Config_Random();
 		Equilibration();
 	}
	
	string t = to_string(temp);
	t.resize(3);
	string field = to_string(h);
	field.resize(3);
  filename = "temp_" + t + "_h_" + field + ".dat";
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

// ************************

void Config_Random(){
	for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
}

// ************************

void Config_Old(){
	ifstream old("Data/config.final");
	for (int i=0; i<nspin; ++i)
  {
  	old >> s[i];
  }
  old.close();
}

// ************************

void Equilibration(){
	for(int i=0; i<100000; ++i)	Move();
	return;
}

// ************************ Move ************************

void Move()
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro == "Metropolis") //Metropolis
    {
			accepted += Metro_Boltzmann(o);
			attempted++;
    }
    if(metro == "Gibbs") //Gibbs sampling
    {
			Gibbs_Boltzmann(o);
    }
  }
}

// ************************ Campionamento ************************

int Metro_Boltzmann(int o){
	
	int new_s = -s[o];																	        	// cambio o-esimo spin
	double Delta_E = Boltzmann(new_s, o) - Boltzmann(s[o], o);		// (energia nuova - energia vecchia)
	double A = min( 1., exp(-beta*Delta_E) );

	double m = rnd.Rannyu();
	if(m<A){																											// se cambio spin, restituisco 1 (---> accepted)
		s[o]=new_s;
		return 1;
	}
	else	return 0;
}

// ************************

void Gibbs_Boltzmann(int o){
	double E_up = Boltzmann(1,o);
	double E_down = Boltzmann(-1,o);
	double Z = exp(-beta*E_up) + exp(-beta*E_down);
	double p = exp(-beta*E_up)/Z;	//probabilità di spin up
	double m = rnd.Rannyu();
	if(m<p)	s[o] = 1;
	else s[o] = -1;
	return;
}

// ************************

double Boltzmann(int sm, int ip)	//calcola energia per il singolo spin (se mi serve la differenza devo chiamarla 2 volte con spin diversi)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;					//somma solo sui primi vicini
  return ene;
}

// ************************ Misura valori istantanei ************************

void Measure()														//valori "istantanei", per energia interna, magnetizzazione etc. mi servono medie!
{
  int bin;
  double u = 0.0, m = 0.0;								// u serve per U e C, m per X e M

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		 m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
 	walker[im] = m;
  walker[ix] = m*m;
}

void Reset(int iblk) //Reset block averages
{
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

// ************************ Accumulo (blocchi) e medie finali ************************

void Accumulate(void) //Update block averages							//NB: in ix e ic non ho ancora valori della X e C, ma solo somme di m^2 e u^2
{
   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

// ************************

void Averages(int iblk) //Print results for current block
{
   ofstream out[4];
   const int wd=12;
      	  
    //cout << "Block number " << iblk << endl;
    out[iu].open("Data/" + metro + "/energy/" + filename,ios::app);
    out[ic].open("Data/" + metro + "/heat_capacity/" + filename,ios::app);
    if(h != 0)	out[im].open("Data/" + metro + "/magnetization/" + filename,ios::app);
    if(h == 0)	out[ix].open("Data/" + metro + "/susceptibility/" + filename,ios::app);		//...altrimenti non apro file
    
    double stima[4];																										//*** singolo blocco ***
    stima[iu] = blk_av[iu]/blk_norm/(double)nspin; 												//< Energia	>									
    double stima_u2 = blk_av[ic]/blk_norm/(double)nspin; 									//< Energia^2	>	
    stima[ic] = beta*beta*( stima_u2 - stima[iu]*stima[iu]*nspin );				//Capacita' termica						ATTENZIONE a /nspin!!!!
    stima[im] = blk_av[im]/blk_norm/(double)nspin; 												//< somma di spin >     ----> Magnetizzazione	
    if(h==0)	stima[ix] = beta*blk_av[ix]/blk_norm/(double)nspin; 				//< (somma di spin)^2 > ----> Susciettivita'
    else			stima[ix] = 0;																							// NB: formula valida solo a campo h==0
    
    for(int i=0; i<4; i++){																							//*** globale all'avanzare dei blocchi ***
   		glob_av[i]  += stima[i];										 						
    	glob_av2[i] += stima[i]*stima[i];
    	err[i] = Error(glob_av[i],glob_av2[i],iblk);
    	out[i] << setw(wd) << iblk <<  setw(wd) << stima[i] << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err[i] << endl;
    	out[i].close();
		}

    //cout << "----------------------------" << endl << endl;
}

// ************************ Risultati ************************

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("Data/config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

	if(metro == "Metropolis" && h==0){
	    ofstream out_acc("Data/Metropolis/acceptance_rate.dat", ios::app);
  	  out_acc << temp << " " << (double)accepted/(double)attempted << endl;
    	out_acc.close();
  }

  rnd.SaveSeed();
}

// ************************ Funzioni supplementari ************************

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

// ************************

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
