CC 	= ${CCP}
ERASE = rm -f
CFLAGS = -g -O3 -ffast-math -Wall -c -I${CFIO_INC_DIR} -I ../
LFLAGS = -g -O3 -ffast-math -Wall -I${CFIO_INC_DIR} -L${CFIODIR} -lcfitsio -lm
OBJECTS = e2adc.o random.o readtext.o

all: e2adc

e2adc: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o e2adc
	cp e2adc ../../bin

random.o: ../ancillary/random.cpp ../ancillary/random.h
	$(CC) $(CFLAGS) ../ancillary/random.cpp -o random.o
readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

e2adc.o: e2adc.cpp
	$(CC) $(CFLAGS) e2adc.cpp -o e2adc.o

clean:
	$(ERASE) $(OBJECTS)
