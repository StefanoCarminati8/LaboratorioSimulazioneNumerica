CC = g++
CFLAGS = -Wall -O3 --std=c++11

main2.exe : main2.o random.o rw.o metodoblocchi.o
	$(CC) random.o main2.o rw.o metodoblocchi.o -o main2.exe
main2.o : main2.cpp
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
rw.o : rw.cpp rw.h
	$(CC) -c rw.cpp -o rw.o $(CFLAGS)
metodoblocchi.o : metodoblocchi.cpp metodoblocchi.h
	$(CC) -c metodoblocchi.cpp -o metodoblocchi.o $(CFLAGS)
clean :
	rm *.o main2.exe seed.out
run :
	make && ./main2.exe
