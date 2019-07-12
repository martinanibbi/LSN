CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_02.exe: es_02.o random.o class_walk.o class_walk.h
	$(CC) random.o class_walk.o es_02.o -o es_02.exe
es_02.o: es_02.cpp random.o class_walk.o
	$(CC) -c es_02.cpp -o es_02.o $(CFLAGS)
class_walk.o: class_walk.cpp class_walk.h
	$(CC) -c class_walk.cpp -o class_walk.o $(CFLAGS)

es_01.exe: es_01.o random.o
	$(CC) random.o es_01.o -o es_01.exe
es_01.o: es_01.cpp random.o
	$(CC) -c es_01.cpp -o es_01.o $(CFLAGS)
	
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
