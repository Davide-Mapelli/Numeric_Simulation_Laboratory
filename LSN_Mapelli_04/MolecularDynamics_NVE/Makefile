CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.cpp
	$(CC) MolDyn_NVE.cpp -o MolDyn_NVE.exe

remove output:
	rm output_solid_ekin.dat
	rm output_solid_epot.dat
	rm output_solid_etot.dat
	rm output_solid_temp.dat

# remove average:
	# rm *.out


clean :
	rm *.o MolDyn_NVE.exe
