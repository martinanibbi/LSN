/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 

  Input();             //Inizialization
  int nconf = 1;			 //serve per frames
  
  int N_blocchi = 100;
	int L = nstep/N_blocchi;			//numero di passi per ogni blocco
	
	ofstream kin_out, pot_out, etot_out, temp_out, pres_out;
	kin_out.open(file_path + "ave_kin.out");
	pot_out.open(file_path + "ave_pot.out");
	etot_out.open(file_path + "ave_etot.out");
	temp_out.open(file_path + "ave_temp.out");
	pres_out.open(file_path + "ave_pres.out");
	
	// 0->kin, 1->pot, 2->etot, 3->temp, 4->pres
	double average[5];
	double average2[5];
	double error[5];
		for(int i=0; i<5; i++)	error[i] = 0;	//al primo blocco l'errore statistico è nullo
	double block[5];
	double block2[5];	
	
 	for(int n=0; n<N_blocchi; n++){			//scorrono i blocchi
	  for(int i=0; i<5; i++){						//svuoto medie del singolo blocco
	  	block[i]=0;		
	  	block2[i]=0;
	  }
  
  	for(int istep=1; istep <= L; ++istep){ //eseguo il blocco (con scrittura in frames)
    	 Move();           //Move particles with Verlet algorithm
   		 Measure();				 //Properties measurement
   		 block[0] += stima_kin;
  		 block[1] += stima_pot;
  		 block[2] += stima_etot;
  		 block[3] += stima_temp;
  		 block[4] += stima_pres;
  		 
  		 if(frames){
  		 	 if(istep%10 == 0){  
					ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
	    		nconf += 1;
	  	 	 }  	
  		 }
  	}
  	for(int i=0; i<5; i++){		//medie sul singolo blocco 
  		 block[i] = (double)(block[i]/L);
  		 block2[i] = (block[i])*(block[i]);			
  	}
  	if(n%10 == 0) cout << "Number of blocks: " << n+10 << endl;
  
 	  for(int i=0; i<5; i++){		//media ed errore all'avanzare dei blocchi
  	   average[i] = (double)((average[i])*n + block[i])/(n+1.);
  		 average2[i] = (double)((average2[i])*n + block2[i])/(n+1.);				
  		 if( n!=0 )	error[i] = sqrt( (average2[i] - (average[i])*(average[i]))/n);
 		}
 
 		kin_out << n+1 << " " << average[0] << " " << error[0] << endl;
  	pot_out << n+1 << " " << average[1] << " " << error[1] << endl;
  	etot_out << n+1 << " " << average[2] << " " << error[2] << endl;
  	temp_out << n+1 << " " << average[3] << " " << error[3] << endl;  
  	pres_out << n+1 << " " << average[4] << " " << error[4] << endl;  
 
  }
  ConfFinal();         //Write final configuration to restart
	kin_out.close();
	pot_out.close();
	etot_out.close();
	temp_out.close();
	pres_out.close();
	
	ofstream last_temp(file_path + "temp_prog.out", ios::app);
	last_temp << average[3] << " " << error[3] << endl;
	last_temp.close();

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
  ReadInput >> iprint;

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

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

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

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

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

void Measure(){ //Properties measurement		//stampa valori istantanei		//valori sono per particella!
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;
	
  Epot.open(file_path + "output_epot.dat",ios::app);			
  Ekin.open(file_path + "output_ekin.dat",ios::app);
  Temp.open(file_path + "output_temp.dat",ios::app);
  Etot.open(file_path + "output_etot.dat",ios::app);
	Pres.open(file_path + "output_pres.dat",ios::app);
	
  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  
//Pressure
	stima_pres = 0.0;
	//sommo su tutte le coppie di particelle
	for(int i=0; i < npart-1; i++){				//prima particella: da 0 a npart - 2
		for(int j=npart-1; j > i; j--){			//seconda particella: da npart - 1 a i + 1		
			double d = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]) );
			stima_pres += ( pow((1./d),12) -0.5*pow((1./d),6) );
		}
	}
	stima_pres *= 16./(double)vol;	//ATTENZIONE: volume, temperatura e densità già corrette se real == true!!
	stima_pres += rho*stima_temp;		//densità * temperatura 	

//correzione per fluidi reali  
  if(real){		
  	stima_temp *= e_kb;
  	stima_kin *= epsilon;
  	stima_pot *= epsilon;
  	stima_pres *= epsilon/(sigma*sigma*sigma);	
  }
  stima_etot = stima_kin + stima_pot; //Total energy per particle (valida sia nel caso reale che in LJUnits)

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
	Pres << stima_pres << endl;
	
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Pres.close();

    return;
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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
