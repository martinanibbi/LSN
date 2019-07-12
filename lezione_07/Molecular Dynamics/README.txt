ATTENZIONE:

in input.dat scrivo:
	1) numero per indicare da quale configurazione partire ---> parto da FCC se <=0, da vecchia configurazione se >0
	2) fase: solid / liquid / gas
	3) elemento: Argon / Krypton / LJUnits
	4) frames: y/n
NB: se problema ad aprire input elemento parto da LJUnits, ma potrebbero esserci problemi se parto dalla vecchia configurazione!!!!

ESERCIZI:
sono presenti 3 programmi in python che cambiano opportunamente i parametri di "Data/input.dat"
	1) per esercizi 1, 2 e 3:
		 execute_LJU_restart.py --> eseguo il programma 1 volta partendo dalla configurazione fcc e N volte riscalando la temperatura
		 FASE: liquid
		 ELEMENTO: LJUnits
		 (ad ogni esecuzione cancello i file di output istantanei, altrimenti dati vengono messi in coda...)
		 
	2) per esercizio 4:
		 execute_real_restart.py --> eseguo il programma per ogni fase ed ogni elemento N volte (partendo da FCC e riscalando la temperatura N 	 
		 volte)
		 execute_real.py --> eseguo il programma per ogni fase ed ogni elemento 1 volta partendo dalla vecchia configurazione.
		 
NOTA:
in tutti i python scripts NON viene eseguito frames: per aggiungere questa funzione modificare a mano il file input.dat
	
	
	
	
