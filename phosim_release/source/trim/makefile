CC 	= ${CCP}
ERASE = rm -f
CFLAGS = -g -O3 -ffast-math -Wall -c -I ../ -I ../../
LFLAGS = -g -O3 -ffast-math -Wall -lm -lz
OBJECTS = main.o trim.o readtext.o
OBJECTS2 = test_trim.o trim.o readtext.o
all: trim test_trim
trim: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o ./trim
	cp trim ../../bin/.
test_trim: $(OBJECTS2)
	$(CC) $(OBJECTS2) $(LFLAGS) -o ./test_trim
readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o
trim.o: trim.cpp
	$(CC) $(CFLAGS) trim.cpp -o trim.o
test_trim.o: test_trim.cpp
	$(CC) $(CFLAGS) test_trim.cpp -o test_trim.o
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp -o main.o
clean:
	$(ERASE) $(OBJECTS)
