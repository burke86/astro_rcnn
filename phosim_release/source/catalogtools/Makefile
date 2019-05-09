CC 	= g++

CFIO_INC_DIR=../cfitsio/include/ 
CFIODIR=../cfitsio/lib/

ERASE = rm -f
//CFLAGS = -g  -ffast-math -Wall -c -I${CFIO_INC_DIR} -I ../
//LFLAGS = -g -ffast-math -Wall  -lm -I${CFIO_INC_DIR} -L${CFIODIR} -lcfitsio -I ../ -L ../ 
CFLAGS = -g  -O3 -ffast-math -Wall -c -I${CFIO_INC_DIR} -I ../
LFLAGS = -g -O3 -ffast-math -Wall  -lm -I${CFIO_INC_DIR} -L${CFIODIR} -lcfitsio -I ../ -L ../ 



OBJECTS = random.o readtext.o galaxycomplex.o sphericalcell.o 

all: phosimcatgen

phosimcatgen: phosimcatgen.o $(OBJECTS)
	 $(CC) -o phosimcatgen  phosimcatgen.o $(OBJECTS) $(LFLAGS)

phosimcatgen.o: phosimcatgen.cpp
	$(CC) $(CFLAGS)  phosimcatgen.cpp -o  phosimcatgen.o

random.o: ../ancillary/random.cpp
	$(CC) $(CFLAGS) ../ancillary/random.cpp -o random.o

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

galaxycomplex.o: galaxycomplex.cpp
	$(CC) $(CFLAGS) galaxycomplex.cpp  -o galaxycomplex.o

sphericalcell.o: sphericalcell.cpp
	$(CC) $(CFLAGS) sphericalcell.cpp  -o sphericalcell.o

clean:
	$(ERASE) $(OBJECTS) phosimcatgen.o
