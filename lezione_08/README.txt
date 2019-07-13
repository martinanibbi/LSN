es_01: legge mu e sigma da file, calcolo energia media con data blocking, salvo mu, sigma e ultimo vaore dell'energia in final_energy.dat.
es_02: trova mu e sigma che minimizzano energia leggendo il file final_energy.dat + campiona |psi|^2 con i valori trovati

execute.py:
	- esegue es_01 molte volte con diversi valori di mu e sigma.
	- esegue es_02 per trovare i mu e sigma migliori e campionare |psi|^2
	- esegue es_01 con mu e sigma migliori, con cui trova energia potenziale al variare dei blocchi
	
QMC_1D: PIGS/PIMC (cambiare direttamente il nome del file di input in qmc1d.cpp)
	NOTA: È stato usato il generatore di numeri casuali di RooT
