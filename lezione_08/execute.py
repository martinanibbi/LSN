import os, sys

#ordine:
#	1) posizione iniziale
#	2) delta
#	3) mu
#	4) sigma

os.system("rm Data/final_energy.dat")		#assicurarsi che funzioni, altrimenti es_02.exe diventa inefficiente

x = '0'
delta = '2.7'
mean = []
std_dev = []
for i in range(0,21):					#divido l'intervallo della media e dev std in 21 parti
	mean.append(str(0.75 + i*0.005))
	#mean.append(str(0.8))
	std_dev.append(str(0.60 + i*0.0025))
	#std_dev.append(str(0.62))

# esercizio 8.01	
#	eseguo più volte il programma es_01.exe cambiando i parametri e salvando i risultati finali in final_energy.dat
for mu in mean:
	for sigma in std_dev:
		file_input = open("Data/input.dat", 'w+')
		file_input.write(x + '\n')
		file_input.write(delta + '\n')
		file_input.write(mu + '\n')
		file_input.write(sigma + '\n')
		file_input.close()
		os.system("./es_01.exe")
		
# esercizio 8.02
# trovo l'energia minima e i parametri mu e sigma nel file final_energy.dat 
# campiono la distribuzione di probabilità: |psi(x)|^2
os.system("./es_02.exe")	
file_final = open("Data/final_GS.dat", 'r')	#leggo mu e sigma che minimizzano l'energia e li ripasso a input.dat
mu = file_final.readline()	#attenzione: con comando readline() già incluso '\n'
sigma = file_final.readline()
file_final.close()

file_input = open("Data/input.dat", 'w+')
file_input.write(x + '\n')
file_input.write(delta + '\n')
file_input.write(mu)	
file_input.write(sigma)
file_input.close()

#eseguo es_01.exe con parametri che minimizzano l'energia 
#		(trovo così espressione dell'energia con data blocking in energy.dat e aggiungo in final_energy.dat una riga con mu, sigma e l'energia 
#		 minima)
os.system("./es_01.exe")	




