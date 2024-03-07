PROG = prog
CC = g++
CPPFLAGS = -g -Wall -std=c++11
OBJS = QuadSeivFac.o

$(PROG) : $(OBJS)
	$(CC) -o $(PROG) $(OBJS)

QuadSeivFac.o : QuadSeivFac.cpp
	$(CC) $(CPPFLAGS) -c QuadSeivFac.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make
