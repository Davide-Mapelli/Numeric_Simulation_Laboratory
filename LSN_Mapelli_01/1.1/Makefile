CC = g++
CFLAGS = -Wall -O3 --std=c++11

main1.exe : main1.o random.o
	$(CC) random.o main1.o -o main1.exe
main1.o : main1.cpp
	$(CC) -c main1.cpp -o main1.o $(CFLAGS)
main2.exe : main2.o random.o
	$(CC) random.o main2.o -o main2.exe
main2.o : main2.cpp
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)
main3.exe : main3.o random.o
	$(CC) random.o main3.o -o main3.exe
main3.o : main3.cpp
	$(CC) -c main3.cpp -o main3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o main1.exe main2.exe main3.exe seed.out
