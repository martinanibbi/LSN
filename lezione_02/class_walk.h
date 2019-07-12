using namespace std;

#ifndef __WalkClass__
#define __WalkClass__

class Walk{
	private:
		int n_steps;											//numero di passi compiuti dall'origine
		
	protected:
		double r;													//distanza dall'origine
		virtual void r_value() = 0;				//r può essere cambiato solo dalle classi figlie (dipende se continuo o discreto)
	
	public:
		Walk();
		virtual ~Walk();
		
		double distance();								//restituisce valore della distanza dall'origine (r)
		virtual void step(Random*) = 0;		//NB: rnd inizializzato nel main, devo passarglielo by reference!
};

// ****************

class Walk_discrete : public Walk{
	private:
		int posizione[3];										//0->x, 1->y, 2->z
		virtual void r_value();							//usata dentro step(), assegna nuovo valore ad r
		
	public:
		Walk_discrete();										//posizione iniziale:origine
		virtual ~Walk_discrete();
		
		int GetX()	{ return posizione[0]; };
		int GetY()	{ return posizione[1]; };
		int GetZ()	{ return posizione[2]; };
		
		virtual void step(Random*);					//viene eseguito un passo in direzione e verso casuale 

};

// ****************

class Walk_continous : public Walk{
	private:
		double posizione[3];								//0->x, 1->y, 2->z
		virtual void r_value();							//usata dentro step(), assegna nuovo valore ad r
		
	public:
		Walk_continous();										//posizione iniziale:origine
		virtual ~Walk_continous();
		
		virtual void step(Random*);					//NB: implementazione sempre attraverso coordinate cartesiane!
};

#endif 
