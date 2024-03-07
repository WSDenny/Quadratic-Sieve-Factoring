PROG = prog
CC = mpic++
CPPFLAGS = -g -Wall
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
