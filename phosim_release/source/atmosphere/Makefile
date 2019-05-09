CFLAGS = -g -O3 -ffast-math -Wall -c -I${CFIO_INC_DIR}  -I${FFTW3_INC_DIR} -I ../ -I ../../
LFLAGS = -g -O3 -ffast-math -Wall -L${FFTW3_DIR} -lfftw3 -L${CFIODIR} -lcfitsio -lm
OBJECTS = main.o atmosphere.o random.o readtext.o
CC  = ${CCS}

all: atmosphere

random.o: ../ancillary/random.cpp ../ancillary/random.h
	$(CC) $(CFLAGS) ../ancillary/random.cpp -o random.o

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp -o main.o

atmosphere: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o atmosphere
	cp atmosphere ../../bin

atmosphere.o: atmosphere.cpp atmosphere.h
	$(CC) $(CFLAGS) atmosphere.cpp -o atmosphere.o

test_atmosphere.o: test_atmosphere.cpp
	$(CCS) $(CFLAGS) test_atmosphere.cpp -o test_atmosphere.o

test_atmosphere: test_atmosphere.o
	$(CC) test_atmosphere.o $(LFLAGS) -o test_atmosphere

clean:
	rm -f $(OBJECTS) atmosphere test_atmosphere
