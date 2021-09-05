CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
remove:
	rm output.ene.0
	rm output.chi.0
	rm output.heat.0
	rm output.mag.0
	rm ave.ene.0
	rm ave.chi.0
	rm ave.heat.0
	rm ave.mag.0
	rm output.ene.1
	rm output.chi.1
	rm output.heat.1
	rm output.mag.1
	rm ave.ene.1
	rm ave.chi.1
	rm ave.heat.1
	rm ave.mag.1

clean :
	rm *.o Monte_Carlo_ISING_1D.exe

