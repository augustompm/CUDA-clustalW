install: bali_score bali_score_elm

clean:
	rm *.o

HEADERS = general.h clustalw.h

CC	= cc
CFLAGS  = -c -g -I$(EXPAT_INC)
LFLAGS	= -g -lm -L$(EXPAT_LIB) -lexpat
EXPAT_LIB	= expat-1.95.6/lib
EXPAT_INC	= expat-1.95.6/include


bali_score : readxml.o init.o util.o bali_score.o
	$(CC) -o $@ readxml.o init.o util.o bali_score.o $(LFLAGS)

bali_score_elm : readxml.o init.o util.o bali_score_elm.o
	$(CC) -o $@ readxml.o init.o util.o bali_score_elm.o $(LFLAGS)


.c.o :	$(HEADERS)
	$(CC) $(CFLAGS) $?

