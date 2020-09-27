CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o random.o metodoblocchi.o
	$(CC) random.o MolDyn_NVE.o metodoblocchi.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h random.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
metodoblocchi.o : metodoblocchi.cpp metodoblocchi.h
	$(CC) -c metodoblocchi.cpp -o metodoblocchi.o $(CFLAGS)

clean :
	rm MolDyn_NVE *.o seed.out
run:
	./clean.sh && make && ./MolDyn_NVE.exe
