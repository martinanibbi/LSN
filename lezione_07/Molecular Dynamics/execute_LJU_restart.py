import os, sys

restart = '0'

for i in range (0,10):
	file_input = open("Data/input.dat", 'w')							# per esercizio 1, 2 e 3 uso fase liquida in unità di Lennard Jones
	file_input.write(restart + "\n")											# prima volta restart == 0 (riparte da FCC), poi riparte da vecchia configurazione
	file_input.write('liquid\n')
	file_input.write('LJUnits\n')
	file_input.write('n' + '\n')													#no frames
	file_input.close()
	os.system("rm -f Data/liquid/LJUnits/*.dat")					# ogni volta che parte il programma cancello file di output istantanei (se esistono)
	os.system("rm -f Data/" + ph + "/" + el + "/output.epot.0")
	os.system("rm -f Data/" + ph + "/" + el + "/output.pres.0")
	os.system("rm -f Data/" + ph + "/" + el + "/output.gave.0")
	os.system("rm -f Data/" + ph + "/" + el + "/output.gofr.0")
		
	os.system("./MolDyn_NVE.exe")
	
	restart = '1'
	
