CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT.exe
remove :
	rm *.0 *.1
