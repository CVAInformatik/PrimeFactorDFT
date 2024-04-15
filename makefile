
CC = g++
CFLAGS = -g 
CPPFLAGS =  -O3

%.o  :  %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@


clean:
	rm *.o


PrimeFactorFFT.o : PrimeFactorFFT.cpp

PrimeFactorDFT.o : PrimeFactorDFT.cpp PrimeFactorDFT.h

PrimeFactorFFT :  PrimeFactorFFT.o PrimeFactorDFT.o 



