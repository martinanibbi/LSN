CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_02.exe: es_02.cpp class_Individual.o class_Sim_Ann.o random.o Travelling_Salesman.h
	mpicxx class_Individual.o class_Sim_Ann.o random.o es_02.cpp -o es_02.exe

es_01.exe: es_01.o class_Individual.o class_Sim_Ann.o random.o Travelling_Salesman.h
	$(CC) random.o class_Individual.o class_Sim_Ann.o es_01.o -o es_01.exe $(CFLAGS)

es_01.o: es_01.cpp class_Individual.o Travelling_Salesman.h random.o 
	$(CC) -c es_01.cpp -o es_01.o $(CFLAGS)

prova.exe: prova.cpp
	mpicxx random.o prova.cpp -o prova.exe

class_Individual.o: class_Individual.cpp Travelling_Salesman.h random.o 
	$(CC) -c class_Individual.cpp -o class_Individual.o $(CFLAGS)
	
class_Sim_Ann.o: class_Sim_Ann.cpp Travelling_Salesman.h random.o 
	$(CC) -c class_Sim_Ann.cpp -o class_Sim_Ann.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o *.exe seed.out
