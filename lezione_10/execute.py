import os, sys

conf = ['circle', 'square']

for c in conf:
	os.system("rm -f Data/" + c + "/*.dat")		
	file_input = open("Data/input.dat", 'w+')
	file_input.write(c + '\n')
	file_input.write('0' + '\n')			# crossover
	file_input.write('0' + '\n')			# swap
	file_input.write('0' + '\n')			# global shift
	file_input.write('0' + '\n')			# local shift
	file_input.write('0' + '\n')			# permutation
	file_input.write('0' + '\n')			# inversion
			# potenza nella scelta dei genitori e numero di iterazioni non servono
			
	file_input.close()
	os.system("./es_01.exe")
	os.system("mpiexec -n 4 ./es_02.exe")
