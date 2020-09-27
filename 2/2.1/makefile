CC = g++
CFLAGS = -Wall -O3 --std=c++11

main2.exe : main2.o random.o metodoblocchi.o funzioni.o Integrale.o
	$(CC) random.o main2.o metodoblocchi.o funzioni.o  Integrale.o -o main2.exe
main2.o : main2.cpp
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
metodoblocchi.o : metodoblocchi.cpp metodoblocchi.h
	$(CC) -c metodoblocchi.cpp -o metodoblocchi.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
Integrale.o : Integrale.cpp Integrale.h
	$(CC) -c Integrale.cpp -o Integrale.o $(CFLAGS)

clean :
	rm *.o main2.exe seed.out
run:
	make && ./main2.exe
