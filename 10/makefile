CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o metodoblocchi.o viaggiatore.o
	$(CC) random.o main.o metodoblocchi.o viaggiatore.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
metodoblocchi.o : metodoblocchi.cpp metodoblocchi.h
	$(CC) -c metodoblocchi.cpp -o metodoblocchi.o $(CFLAGS)
viaggiatore.o : viaggiatore.cpp viaggiatore.h
	$(CC) -c viaggiatore.cpp -o viaggiatore.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.out
run:
	make && ./main.exe
