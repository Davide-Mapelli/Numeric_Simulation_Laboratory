CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o city.o GAs.o 
	$(CC) main.o random.o city.o GAs.o -o main.exe 

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

city.o : city.cpp city.h random.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)

GAs.o : GAs.cpp GAs.h city.h random.h
	$(CC) -c GAs.cpp -o Gas.o $(CFLAGS)
	
clean :
	rm *.o main.exe seed.out
