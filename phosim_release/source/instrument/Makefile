CC 	= ${CCS}
ERASE = rm -f
CFLAGS = -g -O3 -ffast-math -pthread -Wall -std=c++0x -c  -I../
LFLAGS = -g -O3 -ffast-math -pthread -Wall   -lm
OBJECTS = instrument.o random.o instrumentfiles.o readtext.o distortion.o
all: instrument
instrument: $(OBJECTS)
	 $(CC) $(OBJECTS) $(LFLAGS) -o instrument
random.o: ../ancillary/random.cpp ../ancillary/random.h
	$(CC) $(CFLAGS) ../ancillary/random.cpp -o random.o
instrument.o:	instrument.cpp
	$(CC) $(CFLAGS) instrument.cpp -o instrument.o
instrumentfiles.o:	instrumentfiles.cpp
	$(CC) $(CFLAGS) instrumentfiles.cpp -o instrumentfiles.o
readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o
distortion.o:	distortion.cpp
	$(CC) $(CFLAGS) distortion.cpp -o distortion.o
clean:
	$(ERASE) $(OBJECTS)
