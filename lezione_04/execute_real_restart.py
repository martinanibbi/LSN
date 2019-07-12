import os, sys

element = {'Argon', 'Krypton'}
phase = {'solid', 'liquid', 'gas'}

for el in element:
	for ph in phase:	
		restart = '0'
		for i in range (0,4):
			file_input = open("Data/input.dat", 'w')					# per esercizio 4 uso fase solida/liquida/gas per Argon/Krypton
			file_input.write(restart + "\n")									# prima volta restart == 0 (riparte da FCC), poi riparte da vecchia	configurazione
			file_input.write(ph + '\n')
			file_input.write(el + '\n')
			file_input.write('n' + '\n')											#no frames
			file_input.close()
			os.system("rm -f Data/" + ph + "/" + el + "/*.dat")
																												# ogni volta che parte il programma cancello file di output istantanei (se esistono)
			os.system("./MolDyn_NVE.exe")
	
			restart = '1'
