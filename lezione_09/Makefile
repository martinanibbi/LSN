CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe: main.o class_Individual.o class_Population.o random.o Travelling_Salesman.h
	$(CC) random.o class_Individual.o class_Population.o main.o -o main.exe $(CFLAGS)

main.o: main.cpp class_Individual.o Travelling_Salesman.h random.o 
	$(CC) -c main.cpp -o main.o $(CFLAGS)

class_Individual.o: class_Individual.cpp Travelling_Salesman.h random.o 
	$(CC) -c class_Individual.cpp -o class_Individual.o $(CFLAGS)
	
class_Population.o: class_Population.cpp Travelling_Salesman.h random.o 
	$(CC) -c class_Population.cpp -o class_Population.o $(CFLAGS)
	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o *.exe seed.out
