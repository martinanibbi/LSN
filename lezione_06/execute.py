import os, sys

sampling = ['Metropolis', 'Gibbs']
field = ['0', '0.02']
temp = ['0.1']
parameters = ['energy', 'heat_capacity', 'magnetization', 'susceptibility']

Tmin = 0.1
for i in range(1,31):
	temp.append(str( Tmin + i*0.1))
	
os.system("rm -f Data/Metropolis/acceptance_rate.dat")	
for s in sampling:
	for p in parameters:
		os.system("rm -f Data/" + s +"/" + p + "/*.dat")	
	for h in field:
		for t in temp:
			file_input = open("Data/input.dat", 'w+')
			file_input.write(t + '\n') 
			file_input.write('50' + '\n')								#nspin
			file_input.write('1' + '\n')								#J --> exchange interaction
			file_input.write(h + '\n')
			file_input.write(s + '\n')
			file_input.write('20' + '\n')								#nblk
			file_input.write('10000' + '\n')						#nstep
			file_input.write('old' + '\n')							#old -->riparte da vecchia configrazione, altrimenti riparte (casuale)
			if t==Tmin:
				file_input.write('y' + '\n')
			else: 
				file_input.write('n' + '\n')							#equilibrazione: se riparto da temperatura vicina non è necessario equilibrare...
			file_input.close()
			os.system("./Monte_Carlo_ISING_1D.exe")
