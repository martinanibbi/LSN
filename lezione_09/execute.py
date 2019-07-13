import os, sys

conf = ['circle', 'square']

for c in conf:
	os.system("rm -f Data/" + c + "/*.dat")
	file_input = open("Data/input.dat", 'w+')
	file_input.write(c + '\n')
	file_input.write('0.50' + '\n')			# crossover
	file_input.write('0.2' + '\n')			# swap
	file_input.write('0.5' + '\n')			# global shift
	file_input.write('0.125' + '\n')			# local shift
	file_input.write('0.125' + '\n')			# permutation
	file_input.write('0.125' + '\n')			# inversion
	file_input.write('3' + '\n')				# potenza nella scelta dei genitori
	file_input.write('300' + '\n')			# numero di iterazioni
	file_input.close()
	os.system("./main.exe")
