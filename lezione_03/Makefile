CC = g++
CFLAGS = -Wall -O3 --std=c++11

es_01.exe: es_01.o random.o class_price.o
	$(CC) random.o class_price.o es_01.o -o es_01.exe
es_01.o: es_01.cpp random.o
	$(CC) -c es_01.cpp -o es_01.o $(CFLAGS)

class_price.o : class_price.cpp class_price.h
	$(CC) -c class_price.cpp -o class_price.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o es_01.exe seed.out
