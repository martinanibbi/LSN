import os, sys

#only equilibration phase for 1s, 2p with Uniform distribution (it will be used as well for Gauss distribution)

#ordine:
#	1,2,3) posizione iniziale
#	4) delta
#	5) stato (1s/2p)
#	6) distribuzione --> UNIFORME
# 7) n ---> eseguo solo fase di equilibrazione

state = ["1s", "2p"]
delta = []
for i in range(1, 21):
	delta.append(i*0.2)			#vario delta da 0.2 a 4	

position = ['0.5', '0.7', '1.2']	#posizione iniziale ragionevole sia per 1s che per 2p, fissata 	
	
for s in state:
	folder = "Data/" + s + "/Uniform/"
	os.system("rm -f " + folder + "equilibration.dat")	
			#il programma es_01.exe aggiunge una riga ad ogni esecuzione, lo svuoto ad ogni ciclo di equilibration.py
	for d in delta:
		file_input = open("Data/input.dat", 'w+')
		file_input.write(position[0] + '\n')	#posizione iniziale
		file_input.write(position[1] + '\n')
		file_input.write(position[2] + '\n')	
		file_input.write(str(d) + '\n')	#delta (su cui ciclo)
		file_input.write(s + '\n')	#stato atomico
		file_input.write('Uniform' + '\n')
		file_input.write('n' + '\n')	#eseguo solo equilibrazione
		file_input.close()
		os.system("./es_01.exe")
		
#controllo i risultati in equilibration.dat sia per 1s che per 2p e scelgo delta con acettazione ~50% 
#cambio "a mano" i parametri in execution.py
# ---> per lo stato 1s: delta = 1.5
# ---> per lo stato 2p: delta = 3.0

