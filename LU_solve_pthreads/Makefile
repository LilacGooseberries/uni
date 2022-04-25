CFLAGS = -Wall -Werror -O3 -pthread

solve: main.o matrix.o
	g++ main.o matrix.o $(CFLAGS) -lm -o solve

main.o: main.cpp matrix.o
	g++ $(CFLAGS) -c main.cpp

matrix.o: matrix.cpp matrix.h
	g++ $(CFLAGS) -c matrix.cpp

clean:
	rm -f *.o solve
