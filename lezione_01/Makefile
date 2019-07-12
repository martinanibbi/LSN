CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_03.exe: es_03.o random.o
	$(CC) random.o es_03.o -o es_03.exe
es_03.o: es_03.cpp random.o
	$(CC) -c es_03.cpp -o es_03.o $(CFLAGS)

es_02.exe: es_02.o random.o
	$(CC) random.o es_02.o -o es_02.exe
es_02.o: es_02.cpp random.o
	$(CC) -c es_02.cpp -o es_02.o $(CFLAGS)

es_01_3.exe: es_01_3.o random.o
	$(CC) random.o es_01_3.o -o es_01_3.exe
es_01_3.o: es_01_3.cpp random.o
	$(CC) -c es_01_3.cpp -o es_01_3.o $(CFLAGS)

es_01_2.exe: es_01_2.o random.o
	$(CC) random.o es_01_2.o -o es_01_2.exe
es_01_2.o: es_01_2.cpp random.o
	$(CC) -c es_01_2.cpp -o es_01_2.o $(CFLAGS)

es_01_1.exe: es_01_1.o random.o
	$(CC) random.o es_01_1.o -o es_01_1.exe 
es_01_1.o: es_01_1.cpp random.o
	$(CC) -c es_01_1.cpp -o es_01_1.o $(CFLAGS)
	
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
