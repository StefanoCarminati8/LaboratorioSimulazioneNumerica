CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o metodoblocchi.o funzioni.o metropolis.o
	$(CC) random.o main.o metodoblocchi.o funzioni.o metropolis.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
metodoblocchi.o : metodoblocchi.cpp metodoblocchi.h
	$(CC) -c metodoblocchi.cpp -o metodoblocchi.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.out


run:
	make && ./main.exe
