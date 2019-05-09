CC 	= ${CCP}
ERASE = rm -f
CFLAGS = -g -Wall -c -I${FFTW3_INC_DIR}
LFLAGS = -g -Wall -lm -I${FFTW3_INC_DIR} -L${FFTW3_DIR} -lfftw3
OBJECTS = unittest.o ../source/raytrace/raytrace.o

all: unittest

unittest: $(OBJECTS)
	 $(CC) $(OBJECTS) $(LFLAGS) -o unittest

unittest.o: unittest.cpp
	$(CC) $(CFLAGS) unittest.cpp -o unittest.o

clean: 
	$(ERASE) $(OBJECTS)
