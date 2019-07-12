using namespace std;

#ifndef __PriceClass__
#define __PriceClass__


// nelle classi figle cambia solo la valutazione di S(T), calcolo del prezzo dell Call/Put options non cambia!

class Price{
	protected:
		double S0, T, K, r, sigma;
	
	public:
		Price(double s0, double t, double k, double R, double Sigma);
		virtual ~Price();
	
		virtual double asset_price(Random* rnd)=0;									//prezzo di mercato
	
		double Call_price(Random* rnd);										//prezzo dell'opzione in base a UNA valutazione del prezzo di mercato al tempo T
		double Put_price(Random*rnd);											// " "
};

class Price_direct : public Price{
	public:
		Price_direct(double s0, double t, double k, double R, double Sigma);
		virtual ~Price_direct();
		
		virtual double asset_price(Random* rnd);									//prezzo di mercato
};


class Price_discrete : public Price{
	public:
		Price_discrete(double s0, double t, double k, double R, double Sigma);
		virtual ~Price_discrete();
		
		virtual double asset_price(Random* rnd);									//prezzo di mercato
};

#endif 
