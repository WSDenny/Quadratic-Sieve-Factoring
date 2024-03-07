PROG = prog
CC = g++
CPPFLAGS = -g -Wall -fopenmp
LINKFLAGS = -fopenmp
OBJS = QuadSeivFac.o

$(PROG) : $(OBJS)
	$(CC) $(LINKFLAGS) -o $(PROG) $(OBJS)

QuadSeivFac.o : QuadSeivFac.cpp
	$(CC) $(CPPFLAGS) -c QuadSeivFac.cpp

clean:
	rm -f core $(PROG) $(OBJS)

rebuild:
	make clean
	make