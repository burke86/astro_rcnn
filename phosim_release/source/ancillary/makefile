CC = ${CCS}
CFLAGS = -g -Wall -c -I${CFIO_INC_DIR} -I${FFTW3_INC_DIR} -I ../ -I ../../
LFLAGS = -g -Wall -L${FFTW3_DIR} -lfftw3 -L${CFIODIR} -lcfitsio -lm -lz

all:
	make test_readtext test_random

readtext.o: readtext.h readtext.cpp
	$(CC) $(CFLAGS) readtext.cpp -o readtext.o

test_readtext.o: test_readtext.cpp readtext.h readtext.cpp
	$(CC) $(CFLAGS) test_readtext.cpp -o test_readtext.o

test_readtext: test_readtext.o readtext.o
	$(CC) $(LFLAGS) test_readtext.o readtext.o -o test_readtext

random.o: random.h random.cpp
	$(CC) $(CFLAGS) random.cpp -o random.o

test_random.o: test_random.cpp random.cpp random.h
	$(CC) $(CFLAGS) test_random.cpp -o test_random.o

test_random: random.o test_random.o
	$(CC) $(LFLAGS) test_random.o random.o -o test_random

clean:
	rm -f $(OBJECTS) test_readtext
