import os, sys

phase = ['solid', 'liquid', 'gas']

for i in phase:
	os.system("rm Data/" + i + "/*.0")
	file_input = open("Data/input.dat", 'w')
	file_input.write(i + '\n')
	file_input.write('y' + '\n')	#energia e pressione istantanea: y/n
	file_input.close()
	os.system("./Monte_Carlo_NVT.exe")





