/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     
#include <iostream>  
#include <iomanip>   
#include <fstream>      
#include <cmath>        
#include "MolDyn_NVE.h"

//NB: Measure/Accumulate/Averages provengono da Monte_Carlo_NVT.cpp!!! (anche main modificato...)

using namespace std;

int main(){ 

  Input();             //Inizialization
  int nconf = 1;			 //serve per frames
	
	for(int i=0; i<10000; ++i)	Move();	//equilibrazione
	
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();   
      Measure();  
      Accumulate(); //Update block averages

      if(istep%10 == 0){
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }

    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
		

  return 0;
}

//*********************************************************************************//

//leggo input da file preliminare e decido se: 
//	1) iniziare da fcc o da vecchia configurazione da passare a Input_conf(int)
//	2) fase (liquida, solida, gassosa) da passare a Input_phase(string)
//	3) elemento (Argon, Krypton, LJUnits) da passare a Input_element(string)			//in Input_element assegno bool real


void Input(void){			
	ifstream PreInput;
	PreInput.open("Data/input.dat");
	
	int conf;
	string fr;
	PreInput >> conf;
	PreInput >> phase;
	PreInput >> element;
	PreInput >> fr;
	
	folder = "Data/" + phase + "/";
	file_path = folder + element + "/";
	string phase_file = folder + "input." + phase;
	string element_file = "Data/input." + element;
	if(fr == "y")	frames = true;
	else					frames = false;
	
	Random_start();
//NB: ordine di esecuzione --> 1) phase	2) conf	3) element
	Input_phase(phase_file);
	Input_conf(conf);		//riscalo temperatura sempre in unità di Lennard-Jones
	Input_element(element_file);
	
	//Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl; 

	//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
 
  n_props = 2; //Number of observables

//measurement of g(r)
  igofr = 2;
  nbins = 100;
  n_props = n_props + nbins;						//aumento di 100 la dimensione del vettore per i 100 bin
  bin_size = (box/2.0)/(double)nbins;		//raggio massimo della circonferenza che non interseca il box: box/2 -> non prendo tutte le particelle
	
	
	return;
}

//*******************************//

void Random_start(){
 int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream seed_in("seed.in");
   string property;
   if (seed_in.is_open()){
      while ( !seed_in.eof() ){
         seed_in >> property;
         if( property == "RANDOMSEED" ){
            seed_in >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      seed_in.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
return;
}

//*******************************//

void Input_phase(string file_input){ //Prepare all stuff for the simulation
  ifstream ReadInput;
 
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open(file_input.c_str()); //Read input

	cout << "Phase = " << phase << endl;
	
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblk;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();
}

//*******************************//

// nel file Data/input.dat ho letto conf: 
//		se (conf > 0):                                       leggo da Data/(phase)/old.final
//		se (conf <= 0 v problemi in apertura conifg.final) : leggo da Data/config.0 

void Input_conf(int conf){	
	if(conf > 0){
		Input_conf_old();
		return;	
	}	
	else{
		Input_conf_FCC();
		return;
	}	
}

//*******************************//

void Input_conf_FCC(void){

ifstream ReadConf;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("Data/config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

//*******************************//

void Input_conf_old(void){

ifstream ReadFConf, ReadOConf;			//final & old configuration

//Read initial configuration
	string file_fconf = file_path + "old.final";
	string file_oconf = file_path + "old.0";
  ReadFConf.open(file_fconf);
  ReadOConf.open(file_oconf);
  if(ReadOConf.fail() || ReadFConf.fail()){
  	cout << "Error: final configuration not available" << endl;
  	cout << endl << "----------------------------------------" << endl << endl;
  	Input_conf_FCC();
  	return;
  }
  
  cout << "Read initial configuration from file old.final " << endl << endl;
  
  for (int i=0; i<npart; ++i){	//leggo dati final in x,y,z; leggo dati old in xold, yold, zold;
    ReadFConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    
    ReadOConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadFConf.close();
  ReadOConf.close();

//Prepare initial velocities
  cout << "Prepare velocities from old configuration " << endl << endl;
		
	Move();		
   
//ATTENZIONE: in Move r(t): r->rold e genero nuovo r(t+dt) in r, mentre r(t-dt) viene perso dopo aver calcolato le velocità: 	
//		r(t+dt) serve per riscalare velocità e ricalcolare r(t-dt).
//		infine riporto r(t): rold->r 
		
//Riscalo velocità e posizioni vecchie
	
	double sumv2 = 0.0; //velocità quadratica media
  for (int i=0; i<npart; ++i){
  	sumv2 += (vx[i])*(vx[i]) + (vy[i])*(vy[i]) + (vz[i])*(vz[i]);
	}
  sumv2 /= (double)npart;	

  double fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor //NB: temp in unità di LJ
  for (int i=0; i<npart; ++i){
  	vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

		double x_t = xold[i];	 //variabili d'appoggio
		double y_t = yold[i];
		double z_t = zold[i];
		
		xold[i] = Pbc(x[i] - (vx[i])*2.*delta);	//rold corrette con la giusta velocità
		yold[i] = Pbc(y[i] - (vy[i])*2.*delta);
		zold[i] = Pbc(z[i] - (vz[i])*2.*delta);
		 
		x[i] = x_t;
		y[i] = y_t;
		z[i] = z_t;
  }
  return;
}

//*******************************//

//leggo file dell'elemento, se non si apre uso unità di Lennard-Jones e salvo real = false (serve per output)

void Input_element(string file_element){

	ifstream ReadElem;
	ReadElem.open(file_element.c_str());
	if(ReadElem.fail()){	//se problemi ad aprire file parto da LJUnits
		cout << "The program uses Lennard-Jones units " << endl;
		real = false;
		file_path = folder + "LJUnits/";
		return;
	}
		
	ReadElem >> sigma;
	ReadElem >> e_kb;
	ReadElem >> m;
	
	sigma *= nm;			//NB: nel file di input misurato in nanometri!!!
	m *= amu;					//NB: nel file di input misuro in amu!!!
	epsilon = e_kb*kb;
	
	/*
		rho /= (sigma*sigma*sigma);			//faccio tutti i conti in unità di LJ e modifico direttamente l'output in Measure()
		vol *= sigma*sigma*sigma;
		temp *= e_kb;
	*/
	
	cout << "Element = " << element << endl;
	cout << endl << "The program uses IS units " << endl;
 	cout << "Mass of particles = " << m << " kg" << endl;
 	cout << "Energy unity = " << epsilon << " J" << endl;
 	cout << "Space unity = " << sigma << " m" << endl;
 	cout << "Temperature = " << temp*e_kb << " K" << endl;	
	real = true;	
	
	ReadElem.close();
}

//*********************************************************************************//

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

//*********************************************************************************//

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

//*********************************************************************************//

void Measure()
{
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {

// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)            <-------------------- !!!!
		
		 //double N_effective = 0.0;
		 if(dr < box*0.5)															//NB: particelle potrebbero essere più distanti di rmax (devo scartarle anche nel 
		 {																										//conteggio della normalizzazione!)
		 		int indice = int(dr/bin_size) + 2;					//bin: ((ri-rj)/rmax) * nbins + 2
		 		walker[indice] += 2;
		 		//N_effective += 1.;
		 }                                                 

     if(dr < rcut)
     {
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

// contribution to energy and virial
       v += vij;
       w += wij;
     }
    }          
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;
}

//****************************************************************//

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
}

//****************************************************************//

void Accumulate(void) //Update block averages		//NB: anche di g(r)!!!!!
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

//****************************************************************//

void Averages(int iblk) //Print results for current block
{
    
   double r, gdir;
   ofstream Gofr, Gave, Epot, Pres, g_r;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    
    Epot.open(file_path + "output.epot.0",ios::app);
    Pres.open(file_path + "output.pres.0",ios::app);
    Gofr.open(file_path + "output.gofr.0",ios::app);	//media in ogni blocco di g(r)
    Gave.open(file_path + "output.gave.0",ios::app); //media finale di g(r)
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    if(real)	stima_pot *= epsilon;    
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    if(real)	stima_pres *= epsilon/(sigma*sigma*sigma);
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

//Potential energy per particle
    Epot << setw(wd) << iblk <<  " " << setw(wd) << stima_pot << " " << setw(wd) << glob_av[iv]/(double)iblk << " " << setw(wd) << err_pot << endl;
//Pressure
    Pres << setw(wd) << iblk << " " << setw(wd) << stima_pres << " " << setw(wd) << glob_av[iw]/(double)iblk << " " << setw(wd) << err_press << endl;

//g(r)
		
		for(int k=igofr; k<igofr+nbins; ++k)
		{
				r = (k-2.)*bin_size;
				double dvol = (4.*pi/3.) * ( (r+bin_size)*(r+bin_size)*(r+bin_size) - r*r*r );		//serve per normalizzazione (dipende da r)
				gdir = blk_av[k]/(double)(rho*dvol*blk_norm)/(double)npart; 
				glob_av[k] += gdir;
				glob_av2[k] += gdir*gdir;		
		}
		
		if(iblk == nblk)	//stampo solamente all'ultimo...
		{
				for(int k=igofr; k<igofr+nbins; ++k)
				{
						r = (k-2.)*bin_size;
						err_gdir = Error(glob_av[k], glob_av2[k], iblk);
						Gave << r << " " << glob_av[k]/(double)iblk << " " << err_gdir << endl;	
				}
		}
	
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
}


//*********************************************************************************//

void ConfFinal(void){ //Write final and old configuration
  ofstream WriteConf, WriteOld;

  cout << "Print final configuration to file " << file_path << "old.final " << endl << endl;
  WriteConf.open(file_path + "old.final");
  cout << "Print old configuration to file " << file_path << "old.0 " << endl << endl;
  WriteOld.open(file_path + "old.0");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  WriteOld.close();
  return;
}

//*********************************************************************************//

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

//*********************************************************************************//

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

//****************************************************************//

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
