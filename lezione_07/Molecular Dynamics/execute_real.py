import os, sys

element = {'Argon', 'Krypton'}
phase = {'solid', 'liquid', 'gas'}

for el in element:
	for ph in phase:	
		file_input = open("Data/input.dat", 'w')					# per esercizio 4 uso fase solida/liquida/gas per Argon/Krypton
		file_input.write('1' + "\n")												# riparto da vecchia	configurazione
		file_input.write(ph + '\n')
		file_input.write(el + '\n')
		file_input.write('n' + '\n')											#no frames
		file_input.close()
		os.system("rm -f Data/" + ph + "/" + el + "/*.dat")		# ogni volta che parte il programma cancello file di output istantanei (se esistono)
		os.system("rm -f Data/" + ph + "/" + el + "/output.epot.0")
		os.system("rm -f Data/" + ph + "/" + el + "/output.pres.0")
		os.system("rm -f Data/" + ph + "/" + el + "/output.gave.0")
		os.system("rm -f Data/" + ph + "/" + el + "/output.gofr.0")
		
		os.system("./MolDyn_NVE.exe")
	

