CFLAGS = -g -Wall -Werror -O3

all: avl_db

avl_db: main.o avl.o student.o
	g++ main.o avl.o student.o $(CFLAGS) -o avl_db -lm
main.o: main.cpp student.o avl.o
	g++ $(CFLAGS) -c main.cpp

student.o: student.cpp student.h
	g++ $(CFLAGS) -c student.cpp

avl.o: avl.cpp avl.h student.o
	g++ $(CFLAGS) -c avl.cpp

clean:
	rm -f *.o avl_db
