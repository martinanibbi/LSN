import os, sys

#only equilibration phase for 1s, 2p with Uniform distribution (it will be used as well for Gauss distribution)

#ordine:
#	1,2,3) posizione iniziale
#	4) delta
#	5) stato (1s/2p/3d/4f)
#	6) distribuzione (Uniform/Gauss) ---> NB: 3d/4f solo uniforme!!!
# 7) y ---> eseguo tutto il codice 

position = ['0.5', '0.7', '1.2']	#posizione iniziale ragionevole sia per tutti gli stati
state = ["1s", "2p", "3d", "4f"]
delta = ['1.5', '3', '4', '5']	
	#delta per stati 1s, 2p, 3d, 4f (per 1s e 2p possono essere scelti eseguendo equilibration.py e confrontando i risultati in 	
	#equilibration.dat, 3d e 4f scelti "a mano")

#Uniform -> 1s, 2p, 3d, 4f

for i in range(0,4):
	file_input = open("Data/input.dat", 'w+')
	file_input.write(position[0] + '\n')	#posizione iniziale
	file_input.write(position[1] + '\n')
	file_input.write(position[2] + '\n')	
	file_input.write(delta[i] + '\n')			#delta (diversa per ogni stato)
	file_input.write(state[i] + '\n')			#stato atomico
	file_input.write('Uniform' + '\n')
	file_input.write('y' + '\n')					#eseguo tutto il codice
	file_input.close()
	os.system("./es_01.exe")


#Gauss -> 1s, 2p
for i in range(0,2):
	file_input = open("Data/input.dat", 'w+')
	file_input.write(position[0] + '\n')	#posizione iniziale
	file_input.write(position[1] + '\n')
	file_input.write(position[2] + '\n')	
	file_input.write(delta[i] + '\n')			#delta (diversa per ogni stato)
	file_input.write(state[i] + '\n')			#stato atomico
	file_input.write('Gauss' + '\n')
	file_input.write('y' + '\n')					#eseguo tutto il codice
	file_input.close()
	os.system("./es_01.exe")

