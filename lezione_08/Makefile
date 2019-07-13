CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_01.exe :es_01.o random.o
	$(CC) random.o es_01.o -o es_01.exe
es_01.o : es_01.cpp
	$(CC) -c es_01.cpp -o es_01.o $(CFLAGS)
es_02.exe :es_02.o random.o
	$(CC) random.o es_02.o -o es_02.exe
es_02.o : es_02.cpp
	$(CC) -c es_02.cpp -o es_02.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe
