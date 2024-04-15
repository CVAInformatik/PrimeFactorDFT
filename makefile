
CC = g++
CFLAGS = -g 
CPPFLAGS =  -O3

%.o  :  %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@


clean:
	rm *.o

CalcUtil.o : CalcUtil.cpp  CalcUtil.h

PrimeFactorFFT.o : PrimeFactorFFT.cpp

PrimeFactorDFT.o : PrimeFactorDFT.cpp

SlowFFT.o  : SlowFFT.cpp

Calculator.o : Calculator.cpp 
Calulator2E30.o: Calculator2E30.cpp

PrimeFactorFFT :  PrimeFactorFFT.o PrimeFactorDFT.o Calculator.o Calculator3E30.o  SlowFFT.o CalcUtil.o

GrayCodeGenerator.o  : GrayCodeGenerator.cpp  GrayCodeGenerator.h

GrayCodeGenerator : GrayCodeGenerator.o

SHR.o	: SHR.h SHR.cpp

SHR  : SHR.o

