CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o posizione.o
	$(CC) random.o posizione.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
