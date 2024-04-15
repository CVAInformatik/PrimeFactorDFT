/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/
// for performance measurement

#include "CalculatorType.h"
#include "PrimeFactorDFT.h"

#ifdef PERF
#include "windows.h"
#include "profileapi.h"
#endif

#include "PrimeFactorDFT.h"
#include "SlowFFT.h"
#include <iostream>
#include <random>
#include "Calculator.h"
#include "PrimeTable.h"
#include "CalcUtil.h"
#include "Calculator2E30.h"
#include "EllipticCurve.h"


static std::random_device rd;
static std::mt19937 mt(rd());

#define LIMIT 1000000
#define WCOUNT 100


void ClearData(s64 Length, Data* dreal, Data* dimag)
{
    for (s64 i = 0; i < Length; i++) {
        dreal[i] = 0;
        dimag[i] = 0;
    }
}

void test1()
{
    PrimeFactorDFT pf;

    factorSeq  factors;

    std::cout << "Test1 begin " << std::endl;

    //factors.push_back(2);
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    //factors.push_back(11);
    //factors.push_back(19);

    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* real = new Data[pf.Status()];
        Data* imag = new Data[pf.Status()];
        for (int i = 0; i < pf.Status(); i++)
        {
            if (i % 4 == 1) real[i] = 1.0;
            else if (i % 2 == 0) real[i] = 0.0;
            else real[i] = -1.0;
            imag[i] = 0.0;
        }

        pf.forwardFFT(real, imag);
        for (int i = 0; i < pf.Status(); i++)
            std::cout << i << " : " << real[i] << "  " << imag[i] << std::endl;

        pf.InverseFFT(real, imag);
        for (int i = 0; i < pf.Status(); i++)
            std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] real;
        delete[] imag;
    }
    std::cout << "Test1 end " << std::endl << std::endl;

}

/* 'balanced' representation */
#define AT 3.0
#define TA (-3.0)
#define GC 1.0
#define CG (-1.0)

void InitDNA(s64 Length, Data *DNAreal, Data *DNAimag)
{
    s64  len = Length / 8;

    for (s64 i = 0; i < Length; i++) { 
        DNAreal[i] = 0; 
        DNAimag[i] = 0;
    }

    for (s64 i = 0; i < len; i++)    {
        DNAreal[i] = AT;
        DNAreal[i + len] = TA;
        DNAreal[i + len + len] = GC;
        DNAreal[i + len + len + len] = CG;
    }

    // DNAreal now contains #len AT, #len TA, #len GC, #len CG entries,
    // followed by slightly more than 4xlen zeroes
    // we permute the first 4 x len entries, to generates some pseudo-DNA

#define PERMCOUNT 10

    static std::uniform_int_distribution<uint>* dist;
    dist = new std::uniform_int_distribution<uint>(0, (4 * ((int)len)) - 1);
    for (uint x = 0; x < PERMCOUNT; x++)
       for (uint ix = 0; ix < 4 * len; ix++)
       {
           uint i1 = 0;
           i1 = dist->operator()(mt);
           Data t = DNAreal[ix];
           DNAreal[ix] = DNAreal[i1];
           DNAreal[i1] = t;
       };

    

}

uint InitSubDNA(s64 substringLength, s64 Length, Data* subDNAreal, Data* subDNAimag, Data *DNAreal)
{

    for (s64 i = 0; i < Length; i++) {
        subDNAreal[i] = 0;
        subDNAimag[i] = 0;
    }

    static std::uniform_int_distribution<uint>* dist;

    dist = new std::uniform_int_distribution<uint>(0, (4 * ((int) (Length/8))) - 1);
    uint i1 = 0;
    do {
        i1 = dist->operator()(mt);
    } while (i1 >= ((4 * ((int)(Length / 8))) - substringLength));
    std::cout << "start:  " << i1 << std::endl;
    /* we swap direction ! */
    s64 iy = Length - 1;
    for (uint ix = 0; ix < substringLength; ix++)
        subDNAreal[iy--] = DNAreal[ix + i1];

    return i1;

}



void Multiply(s64 Length , Data *AxBreal, Data* AxBimag, Data* Areal, Data* Aimag, Data* Breal, Data*  Bimag)
{
    for (s64 i = 0; i < Length; i++) {
        Data tr, ti; 

        tr = Areal[i] * Breal[i] - Aimag[i] * Bimag[i];
        ti = Breal[i] * Aimag[i] + Areal[i] * Bimag[i];
        AxBreal[i] = tr;
        AxBimag[i] = ti;
    }

}

void test2DNA()
{

#define SUBSTRINGLENGTH 300

    PrimeFactorDFT pf;

    factorSeq  factors;


    std::cout << "Test2DNA begin " << std::endl;
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    factors.push_back(11);
    factors.push_back(19);
    factors.push_back(31);

    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* DNAreal    = new Data[pf.Status()];
        Data* DNAimag    = new Data[pf.Status()];
        Data* subDNAreal = new Data[pf.Status()];
        Data* subDNAimag = new Data[pf.Status()];
        Data* Matchreal  = new Data[pf.Status()];
        Data* Matchimag  = new Data[pf.Status()];

        InitDNA(pf.Status(),DNAreal, DNAimag);
        InitSubDNA(SUBSTRINGLENGTH, pf.Status(), subDNAreal, subDNAimag, DNAreal );


        pf.forwardFFT(DNAreal, DNAimag);
        pf.forwardFFT(subDNAreal, subDNAimag);

        Multiply(pf.Status(), Matchreal, Matchimag, DNAreal, DNAimag, subDNAreal, subDNAimag);
        pf.ScaledInverseFFT(Matchreal, Matchimag);

        Data max = 0;
        Data max2 = 0;
        s64 maxIndex = 0;
        s64 maxIndex2 = 0;
        for (int i = 0; i < pf.Status(); i++) {
            Data val = (Matchreal[i] * Matchreal[i]) + (Matchimag[i] * Matchimag[i]);
            //std::cout << i << " : " << val << std::endl;
            if (val > max) {
                max2 = max;
                max = val;
                maxIndex2 = maxIndex;
                maxIndex = i;
            }
        }
        std::cout << " maxIndex  : " << maxIndex   << " val : " << max  << std::endl;
        std::cout << " maxIndex2 : " << maxIndex2  << " val : " << max2 << std::endl;
        //pf.InverseFFT(real, imag);
        //for (int i = 0; i < pf.Status(); i++)
        //    std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] DNAreal;
        delete[] DNAimag;
        delete[] subDNAreal;
        delete[] subDNAimag;
        delete[] Matchreal;
        delete[] Matchimag;

    }
    std::cout << "Test2DNA End " << std::endl << std::endl;

}

void test3Convolution()
{


    PrimeFactorDFT pf;

    factorSeq  factors;


    std::cout << "Test3Convolution begin " << std::endl;
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    factors.push_back(11);
    //factors.push_back(19);
    //factors.push_back(31);

    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* A1real = new Data[pf.Status()];
        Data* A1imag = new Data[pf.Status()];
        Data* A2real = new Data[pf.Status()];
        Data* A2imag = new Data[pf.Status()];
        Data* Resreal = new Data[pf.Status()];
        Data* Resimag = new Data[pf.Status()];

        ClearData(pf.Status(), A1real, A1imag);
        ClearData(pf.Status(), A2real, A2imag);
        ClearData(pf.Status(), Resreal, Resimag);

        A1real[0] = 1.0/15;
        A1real[1] = 12.0/15;
        A1real[2] = 1.0/15;
        A1real[3] = 0.5/15;
        A1real[pf.Status() - 1] = 0.5/15;
        A2real[500] = 1.0;

        pf.forwardFFT(A1real, A1imag);
        pf.forwardFFT(A2real, A2imag);
#define COUNT 5
        Multiply(pf.Status(), Resreal, Resimag, A1real, A1imag, A2real, A2imag);

        for (int i = 1; i < COUNT; i++)
            Multiply(pf.Status(), Resreal, Resimag, A1real, A1imag, Resreal, Resimag);
        
        pf.ScaledInverseFFT(Resreal, Resimag);

        //pf.InverseFFT(Resreal, Resimag);
        for (int i = 0; i < pf.Status(); i++)
            std::cout << i << " " << Resreal[i] << "  " << Resimag[i] << std::endl;

        delete[] A1real ;
        delete[] A1imag ;
        delete[] A2real ;
        delete[] A2imag ;
        delete[] Resreal;
        delete[] Resimag;

    }
    std::cout << "Test3Convolution End " << std::endl << std::endl;

}


void test4()
{
    PrimeFactorDFT pf;
    SlowFFT sft;

    factorSeq  factors;

    std::cout << "Test1 begin " << std::endl;

    /*factors.push_back(2);
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);*/
    factors.push_back(11);
    factors.push_back(13);
    factors.push_back(19);
    factors.push_back(31);
    pf.SetFactors(factors);
    sft.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* real = new Data[pf.Status()];
        Data* imag = new Data[pf.Status()];
        Data* sreal = new Data[pf.Status()];
        Data* simag = new Data[pf.Status()];
        for (int i = 0; i < pf.Status(); i++)
        {
            if (i % 4 == 1) sreal[i] = real[i] = 1.0;
            else if (i % 2 == 0) sreal[i] = real[i] = 0.0;
            else sreal[i] = real[i] = -1.0;
            simag[i] = imag[i] = 0.0;
        }
        std::cout << "PFA begin" << std::endl;
        pf.forwardFFT(real, imag);
        std::cout << "PFA end" << std::endl;

        std::cout << "SFT begin" << std::endl;
        sft.forwardFFT(sreal, simag);
        std::cout << "SFT end" << std::endl;

        for (int i = 0; i < pf.Status(); i++) {
            std::cout << "PFA :" << i << " : " << real[i] << "  " << imag[i] << std::endl;
            std::cout << "Ref :" << i << " : " << sreal[i] << "  " << simag[i] << std::endl << std::endl;
        }
        //pf.InverseFFT(real, imag);
        //for (int i = 0; i < pf.Status(); i++)
        //    std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] real;
        delete[] imag;
        delete[] sreal;
        delete[] simag;
    }
    std::cout << "Test1 end " << std::endl << std::endl;

}

void test5perf()
{
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif

    PrimeFactorDFT pf;

    factorSeq  factors;

    std::cout << "Test5 begin " << std::endl;

    //factors.push_back(2);
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    factors.push_back(11);
    //factors.push_back(13);
    //factors.push_back(17);
    factors.push_back(19);
    factors.push_back(31);
    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* real = new Data[pf.Status()];
        Data* imag = new Data[pf.Status()];
        for (int i = 0; i < pf.Status(); i++)
        {
            if (i % 4 == 1)  real[i] = 1.0;
            else if (i % 2 == 0)  real[i] = 0.0;
            else real[i] = -1.0;
            imag[i] = 0.0;
        }
        std::cout << "PFA begin" << std::endl;
#ifdef PERF
        QueryPerformanceFrequency(&Frequency);
        QueryPerformanceCounter(&StartingTime);
#endif
        pf.forwardFFT(real, imag);
#ifdef PERF
        QueryPerformanceCounter(&EndingTime);
        ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
        ElapsedMicroseconds.QuadPart *= 1000000;
        ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
        std::cout << "Length " << pf.Status() << "  Elapsed time (microseconds): " << ElapsedMicroseconds.QuadPart << std::endl;
#endif
        std::cout << "PFA end" << std::endl;


       /* for (int i = 0; i < pf.Status(); i++) {
            std::cout << "PFA :" << i << " : " << real[i] << "  " << imag[i] << std::endl;
        }*/
        //pf.InverseFFT(real, imag);
        //for (int i = 0; i < pf.Status(); i++)
        //    std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] real;
        delete[] imag;
    }
    std::cout << "Test5 end " << std::endl << std::endl;

}

#define PRINT( c) do { std::cout << *c.ItoA() << std::endl;} while (0)


void test6Calc()
{
    CALCULATOR c;
    //char num[] = "261261924691694619461924691649164";
    char num[] = "90000000";
    std::cout << num << std::endl;
    c.Push(num);
    c.Dup();
    PRINT(c);

    //char num1[] = "1000000";
    char num1[] = "723972359729357927536526515164159069889121"
        //"72397235972935792753652651516415906988912444" 
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        ;

    std::cout << num1 << std::endl;
    c.Push(num1);
    c.Dup();
    PRINT(c);
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
    c.Dup();
    c.Mul();
    c.Dup();
    PRINT(c);
}



void test7Calc()
{
    CALCULATOR c;
    char num[] = "261261924691694619461924691649164";
    c.Push(num);
    c.PopStore("test");
    PRINT(c);
    c.PushStore("test");
    PRINT(c);
    c.PushStore("test1");
    c.Push(0xFFFFFFFF);
    PRINT(c);
    c.Push(0x7FFFFFFF);
    PRINT(c);
}

void test8Calc()
{
    CALCULATOR c;

    char num1[] = "500000000";
    char num2[] = "500000001";
    char num3[] = "3996879";
    char num4[] = "4637923";
    char num3000[] = "3996879000";
    char num4000[] = "4637923000";
    char num[] = "261261924691694619461924691649164";


    c.Push(num1);
    c.Push(num2);
    c.Add();
    PRINT(c);
    c.Push(num1);
    c.ChangeSign();
    c.Push(num2);
    c.Add();
    PRINT(c);
    c.Push(num2);
    c.ChangeSign();
    c.Push(num1);
    c.Add();
    PRINT(c);

    c.Push(num3);
    c.Push(num4);
    c.ChangeSign();
    c.Add();
    PRINT(c);

    c.Push(num4);
    c.Push(num3);
    c.ChangeSign();
    c.Add();
    PRINT(c);

    c.Push(num3000);
    c.Push(num4000);
    c.ChangeSign();
    c.Add();
    PRINT(c);

    c.Push(num4000);
    c.Push(num3000);
    c.ChangeSign();
    c.Add();
    PRINT(c);



    c.Push(num);
    c.Dup();
    c.Dup();
    PRINT(c);
    c.Add();
    PRINT(c);
    c.Push(num);
    c.Dup();
    c.Dup();
    PRINT(c);
    c.ChangeSign();
    c.Dup();
    PRINT(c);
    c.Add();
    PRINT(c);
    c.Push(num);
    c.Push(-1);
    c.Add();
    PRINT(c);
    c.Push(1);
    c.Push(num);
    c.ChangeSign();
    c.Add();
    PRINT(c);
    std::cout  << std::endl;
}

void test9Calc()
{
    CALCULATOR c;
    char num1[] = "77777777777777777777777777777777";
    char num2[] = "33333333333333333333333333333333";
    c.Push(num1);
    c.Push(num2);
    c.Add();
    PRINT(c);
    c.Push(num1);
    c.Push(num2);
    c.ChangeSign();
    c.Add();
    PRINT(c);
    c.Push(num2);
    c.Push(num1);
    c.ChangeSign();
    c.Add();
    PRINT(c);
    c.Push(num2);
    c.ChangeSign();
    c.Push(num1);
    c.ChangeSign();
    c.Add();
    PRINT(c);
    c.Push(num2);
    c.ChangeSign();
    c.Push(num1);
    c.Add();
    PRINT(c);
    c.Push(num1);
    c.ChangeSign();
    c.Push(num2);
    c.Add();
    PRINT(c);
}

void test10Calc()
{
    //CALCULATOR c;
    //c.Push(1);
    //c.Push(1);
    //for (int i = 0; i < 1000; i++){
    //    c.Add();
    //    c.Dup();
    //    c.Dup();
    //    std::string* s = c.ItoA();
    //    std::cout << "size " << s->length() << std::endl;
    //    std::cout << *s << std::endl << std::endl;
    //}
    //std::cout  << std::endl << std::endl;

    CALCULATOR c;
    c.Push(1);
    c.Push(1);
    c.Mul();
    PRINT(c);
    c.Push(0);
    c.Push(1);
    c.Mul();
    PRINT(c);
    c.Push(2);
    for (int i = 0; i < 20; i++) {
        c.Dup();
        c.Mul();
        c.Dup();
        std::string* s = c.ItoA();
        std::cout << "size " << s->length() << std::endl;
        std::cout << *s << std::endl << std::endl;
    }
/*    */
}

void test11CalcGCD()
{
    CALCULATOR c;
//    char num[] = "693";//"555555555";
//    char num1[] = "609";//"555";
    char num1[] = "555555555";
    char num[] = "555";
    c.Push(num1);
    c.Push(num);
    c.GCD();
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    char num2[] = "5555";
    char num3[] = "557";
    c.Push(num2);
    c.Push(num3);
    c.GCD();
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    char num4[] = "110893008";
    char num5[] = "7448755608";
    c.Push(num4);
    c.Push(num5);
    c.GCD();
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;

}

void test12CalcSmall()
{
    CALCULATOR c;
    //    char num[] = "693";//"555555555";
    //    char num1[] = "609";//"555";
    char num1[] = "777";
    char num[] = "555";
    c.Push(num1);
    c.Push(num);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num01[] = "77777777777777";
    char num0[] = "555";
    c.Push(num01);
    c.Push(num0);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num21[] = "777";
    char num2[] = "5555555555555555";
    c.Push(num21);
    c.Push(num2);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num31[] = "77777777777777777777777777777777777777";
    char num3[] = "555555555555555555555555555555555555555";
    c.Push(num31);
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num3);
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num41[] = "7777777777777777777777777777777777777777777777777777777777777777777777777777";
    char num4[] = "555555555555555555555555555555555555555555555555555555555555555555555555555555";
    c.Push(num41);
    c.Push(num4);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;

}

void test13QuotientReminder()
{
    CALCULATOR c;
    
    char num1[] = "777777777777777";
    char num[] = "5555555555555";

    //c.Push(num1);
    //c.Push(num);
    //c.QuotientRemainder();
    //std::cout << "Reminder: " <<  *c.ItoA() << std::endl;
    //std::cout << "Quotient: " <<  *c.ItoA() << std::endl;

    //char num21[] = "777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777";
    //char num2[] = "5555555555555555555555555555555555555555555555555555555555555555555555555555555";
    //c.Push(num21);
    //c.Push(num2);
    //c.QuotientRemainder();
    //std::cout << "Reminder: " << *c.ItoA() << std::endl;
    //std::cout << "Quotient: " << *c.ItoA() << std::endl;

#define LIMIT1 100000
#define STEP  177395769

    c.Push(num1);
    c.PopStore("num1");
    c.Push(num);
    c.PopStore("num");
    for (int i = 0; i < LIMIT1; i++)
    {
        c.PushStore("num1");
        //c.Dup();
        //std::cout << "num1: " << *c.ItoA() << std::endl;
        c.PushStore("num");
        c.Push(i);
        c.Push(STEP);
        c.Mul();
        c.Add();
        //c.Dup();
        //std::cout << "num:      " << *c.ItoA() << std::endl;
        c.QuotientRemainder();
        c.PopStore("Remainder:");
        c.PopStore("Quotient: ");
        //c.PushStore("Remainder:");
        //std::cout << "Reminder: " << *c.ItoA() << std::endl;
        //c.PushStore("Quotient: ");
        //std::cout << "Quotient: " << *c.ItoA() << std::endl;
        c.PushStore("num");
        c.Push(i);
        c.Push(STEP);
        c.Mul();
        c.Add();
        c.PushStore("Quotient: ");
        std::cout << "Quotient: " << *c.PrintTOS() << std::endl;
        c.Mul();
        c.PushStore("Remainder:");
        c.Add();
        c.PushStore("num1");
        if (!c.IsEqual()) {
            c.Swap();
            std::cout << "num * Quotient + Remainder: " << *c.ItoA() << std::endl;
            std::cout << "Num1                      : " << *c.ItoA() << std::endl;
        }

    }
    c.ClearStore();
}

void test14Exp()
{
    CALCULATOR c;
    std::string* s;
    c.Push(2);
    c.Push(8624);
    c.Exp();
    c.Push(-1);
    c.Add();
    s = c.ItoA();
    std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;


    c.Push(2);
    c.Push(86243);
    c.Exp();
    c.Push(-1);
    c.Add();
    s = c.ItoA();
    std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;

    c.Push(2);
    c.Push(13466917);
    c.Exp();
    c.Push(-1);
    c.Add();
    s = c.ItoA();
    std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;

    c.Push(2);
    c.Push(20996011);
    c.Exp();
    c.Push(-1);
    c.Add();
    s = c.ItoA();
    std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;

}

void test15Jacobi()
{
	CALCULATOR c;
	//std::cout << "Jacobi ( ";
	//c.Push(9907);
	//c.Dup();
	//std::cout << *c.ItoA() << "/ ";
	//c.Push(1001);
	//c.Dup();
	//std::cout << *c.ItoA() << " ) = ";
	//c.Jacobi();
	//std::cout << *c.ItoA() << std::endl;


	for (int j = 3; j < 60; j += 2)
		for (int i = 1; i < 31; i++) {
			std::cout << "Jacobi ( ";
			c.Push(j);
			c.Dup();
			std::cout << *c.ItoA() << "/ ";
			c.Push(i);
			c.Dup();
			std::cout << *c.ItoA() << " ) = ";
			c.Jacobi();
			std::cout << *c.ItoA() << std::endl;

		}

}


//bool MillerRabin(BInt& number);
//bool MillerRabin(BInt& number, const std::vector<unsigned int>& witnesses);


void test16MillerRabin(char c[]) {

    PrimeTable pt(LIMIT);

    std::vector<unsigned int>  witnesses;

    witnesses.push_back(2);
    unsigned int i = 3;
    while (witnesses.size() < WCOUNT) {
        if (pt.IsPrime(i)) witnesses.push_back(i);
        i = i + 2;
    }

    CALCULATOR cal; 

    /* debug code to be removed*/
    //char YY[] = "4611686014132420609";
    //char m[] =  "2147483649";
//    char YY[] = "2605843007066210305";
//    char m[]   = "1073741824";
    //cal.Push(YY);
    //cal.Push(m);
    //cal.QuotientRemainder();
    //std::cout << "Remainder: " << *cal.ItoA() << std::endl;
    //std::cout << "Quotient : " << *cal.ItoA() << std::endl;
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif
#ifdef PERF
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
#endif


    cal.Push((char*) c );
    for (int itx = 0; itx < 2000; itx++) {
        cal.Push((char*)c);
        cal.Push(2);
        cal.Push(itx);
        cal.Mul();
        cal.Add();    
        BINT temp;
        cal.Pop(temp);
        if (MillerRabin(temp, witnesses)) {
            //cal.Push(temp);
            //std::cout << "probably prime: " << *cal.ItoA() << std::endl;
        }
    }

#ifdef PERF
        QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMicroseconds.QuadPart *= 1000000;
    ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
    std::cout << "Elapsed time(microseconds) : " << ElapsedMicroseconds.QuadPart << std::endl;
#endif
}

void test17()
{
    CALCULATOR c;

    c.Push((char *) "9365865165198658618561865816581658");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }
    c.Push((char*)"999999999");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }
    c.Push((char*)"999");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }

    c.Push((char*)"99");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }

    c.Push((char*)"9");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }
    c.Push((char*)"0");
    c.Dup();
    c.Dup();
    std::cout << "Arg    " << *c.ItoA() << std::endl;
    c.Rand();
    std::cout << "Rand() " << *c.ItoA() << std::endl;

}

void test18()
{
    CALCULATOR c;

    c.Push((char*)"9999999999999999999999999999999999");
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Dup();
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Dup();
    c.ChangeSign();
    c.Swap();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.ChangeSign();
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Push((char*)"999999");
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.ChangeSign();
    c.Push((char*)"999999");
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Push((char*)"999999");
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.ChangeSign();
    c.Push((char*)"999999");
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.Dup();
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.Dup();
    c.ChangeSign();
    c.Swap();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.ChangeSign();
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;

}


void test19()
{
    CALCULATOR cal;
    cal.Push(859433);
    cal.Push(3021377);
    cal.Mul();
    std::string *s = cal.ItoA();
    Factoring((char *) s->c_str());

    Factoring((char*)"19777122841");

}


void test20() 
{
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif

    CALCULATOR cal;
    BINT P;
    BINT A;
    BINT Res;
#define P224 1
#if P224
    // NIST P-224 
#ifdef PERF
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
#endif

    char Pascii[] = "26959946667150639794667015087019630673557916260026308143510066298881";
    cal.Push(Pascii);
    cal.Pop(P);
    cal.Push(2021);
    cal.Dup();
    cal.Dup();
    cal.Mul();
    cal.Mul();
    cal.Push(2021);
    cal.Push(-3);
    cal.Mul();
    cal.Add();
    char Aascii[] = "18958286285566608000408668544493926415504680968679321075787234672564";
    cal.Push(Aascii);
    cal.Add();
    cal.Pop(A);
    SquareRootModM(Res, A, P);
#ifdef PERF
    QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMicroseconds.QuadPart *= 1000000;
    ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
    std::cout << "Elapsed time(microseconds) : " << ElapsedMicroseconds.QuadPart << std::endl;
#endif

    cal.Push(A);
    std::cout << std::endl << " A:   " << *cal.ItoA() << std::endl;
    cal.Push(P);
    std::cout << " P:   " << *cal.ItoA() << std::endl;
    cal.Push(Res);
    cal.Dup();
    std::cout << " Res: " << *cal.ItoA() << std::endl;
    cal.Square();
    cal.Dup();
    std::cout << " Res * Res  : " << *cal.ItoA() << std::endl;
    cal.Push(P);
    cal.Mod();
    std::cout << " Res * Res  mod P: " << *cal.ItoA() << std::endl;
    return;

#else 
    cal.Push((char* )"2147483647");
    //cal.Push((char* )"43");
    cal.Pop(P);
    cal.Push((char* ) "3497491");
    //cal.Push((char* ) "6");
    cal.Pop(A);
    cal.Push(P);
    cal.Push(A);
    cal.Jacobi();
    if (cal.IsEqual(1)) {
        SquareRootModPrime(Res, A, P); 
        cal.Push(A);
        std::cout << " A:   " << *cal.ItoA() << std::endl;
        cal.Push(P);
        std::cout << " P:   " << *cal.ItoA() << std::endl;
        cal.Push(Res);
        cal.Dup();
        std::cout << " Res: " << *cal.ItoA() << std::endl;
        cal.Square();
        cal.Dup();
        std::cout << " Res * Res  : " << *cal.ItoA() << std::endl;
        cal.Push(P);
        cal.Mod();
        std::cout << " Res * Res  mod P: " << *cal.ItoA() << std::endl;
        return;
    }
#endif
}


void test21()
{
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif

    CALCULATOR cal;
    BINT res;

    int arg = 200;

#ifdef PERF
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
#endif

    Faculty(res, arg);
#ifdef PERF
    QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMicroseconds.QuadPart *= 1000000;
    ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
    std::cout <<  "Elapsed time(microseconds) : " << ElapsedMicroseconds.QuadPart << std::endl;
#endif
    cal.Push(res);
    std::cout << "res  " << *cal.ItoA() << std::endl;
}

void test22()
{
    Calculator cal;
    BInt temp;
    BInt2E30 temp2;

    cal.Push(2);
    cal.Push(1000);
    cal.Exp();
    cal.Push(-1);
    cal.Add();
    cal.Dup();
    std::cout << "temp " << *cal.ItoA() << std::endl;
 /*   cal.Pop(temp);
    Convert10E9to2E30(temp2, temp);
    for (int i = 0; i < temp2.number.size(); i++) printf(" %3d :  %08X \n", i, temp2.number[i]);*/

    //Convert2E30to10E9(temp, temp2);
    //cal.Push(temp);
    //std::cout << "temp " << *cal.ItoA() << std::endl << std::endl;

    MersenneBInt2E20(temp2, 86243);
    Convert2E30to10E9(temp, temp2);
    cal.Push(temp);
    std::cout << "temp " << *cal.ItoA() << std::endl  << std::endl;

    MersenneBInt2E20(temp2, 209960);
    Convert2E30to10E9(temp, temp2);
    cal.Push(temp);
    std::cout << "temp " << *cal.ItoA() << std::endl << std::endl;

    MersenneBInt2E20(temp2, 209960 *5);
    Convert2E30to10E9(temp, temp2);
    cal.Push(temp);
    std::cout << "temp " << *cal.ItoA() << std::endl << std::endl;

}

void test23()
{
    CALCULATOR cal;
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif



    cal.Push((char *) "293619494913165487615161600641311312659");
    std::cout << " : " << *cal.ItoA() << std::endl;

    std::string* s;

#define astr "1234567890"
    char cstring[] =
        astr astr astr astr astr \
        astr astr astr astr astr \
        astr astr astr astr astr \
        astr astr astr astr astr;

    char* cptr = cstring;
    while (*cptr) {
        cal.Push(cptr);
        printf("               Exp  %s\n", cptr);
        s = cal.ItoA();
        printf("length: %4d   Exp  %s\n", (int)s->length(), s->c_str());
        cptr++;
    }


    for (int i = 0; i < 20; i++) {
        cal.Push(100001);
        cal.Push(i + 1);
        cal.Exp();
        s = cal.ItoA();
        printf("length:  %4d   Exp  %s\n", (int)s->length(), s->c_str());
        //std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;
    }

    cal.Push(2);
    cal.Push(86);
    cal.Exp();
    cal.Push(-1);
    cal.Add();
    s = cal.ItoA();
    std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;

    for (int i = 0; i < 50; i++) {
        cal.Push(2);
        cal.Push(862);
        cal.Push(i);
        cal.Add();
        cal.Exp();
        //cal.Rand();
        s = cal.ItoA();
        std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;
    }

    cal.Push(2);
    cal.Push(8624);
    cal.Exp();
    cal.Push(-1);
    cal.Add();
#ifdef PERF
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
#endif

    s = cal.ItoA();
#ifdef PERF
    QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMicroseconds.QuadPart *= 1000000;
    ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
    std::cout << "Elapsed time(microseconds) : " << ElapsedMicroseconds.QuadPart << std::endl;
#endif

    std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;


    cal.Push(2);
    cal.Push(86243);
    cal.Exp();
    cal.Push(-1);
    cal.Add();
#ifdef PERF
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
#endif

    s = cal.ItoA();
#ifdef PERF
    QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMicroseconds.QuadPart *= 1000000;
    ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
    std::cout << "Elapsed time(microseconds) : " << ElapsedMicroseconds.QuadPart << std::endl;
#endif
    std::cout << "length:   " << s->length() << "   Exp " << *s << std::endl;


}
/*

   Polynomial multiplication  modulo X^7 -1

*/

void test24()
{
    PrimeFactorDFT pf;

    factorSeq  factors;

    std::cout << "Test24 begin " << std::endl;

    factors.push_back(7);

    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* real = new Data[pf.Status()];
        Data* imag = new Data[pf.Status()];

        Data* real1 = new Data[pf.Status()];
        Data* imag1 = new Data[pf.Status()];

        for (int i = 0; i < pf.Status(); i++) {
            real[i] = real1[i] = imag[i] = imag1[i] = 0.0;
        }

        real[0] = 3;
        real[1] = 4;
        real[2] = 2;
        real[3] = 1;

        real1[0] = 2;
        real1[1] = 5;
        real1[2] = 4;
        real1[3] = 2;
        real1[4] = 1;

        pf.forwardFFT(real, imag);
        pf.forwardFFT(real1, imag1);


        for (int i = 0; i < pf.Status(); i++) {
            double re = (real[i] * real1[i]) - (imag[i] * imag1[i]);
            double im = (real[i] * imag1[i]) + (real1[i] * imag[i]) ;
            real[i] = re;
            imag[i] = im;
        }


        pf.ScaledInverseFFT(real, imag);


        for (int i = 0; i < pf.Status(); i++)
            std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] real;
        delete[] imag;
        delete[] real1;
        delete[] imag1;
    }

    for (int j = 0; j < 7; j++) {
        std::cout << std::endl << "j  " << j << std::endl << std::endl;
        if (pf.Status() > 0) {

            Data* real = new Data[pf.Status()];
            Data* imag = new Data[pf.Status()];

            Data* real1 = new Data[pf.Status()];
            Data* imag1 = new Data[pf.Status()];

            for (int i = 0; i < pf.Status(); i++) {
                real[i] = real1[i] = imag[i] = imag1[i] = 0.0;
            }

            real[0] = 3;
            real[1] = 4;
            real[2] = 2;
            real[3] = 1;

            real1[j % 7] = 1;
            real1[(j + 1) % 7] = 1;
            real1[(j + 2) % 7] = 1;

            pf.forwardFFT(real, imag);
            pf.forwardFFT(real1, imag1);


            for (int i = 0; i < pf.Status(); i++) {
                double re = (real[i] * real1[i]) - (imag[i] * imag1[i]);
                double im = (real[i] * imag1[i]) + (real1[i] * imag[i]);
                real[i] = re;
                imag[i] = im;
            }


            pf.ScaledInverseFFT(real, imag);

            for (int i = 0; i < pf.Status(); i++)
                std::cout << i << " " << real[i] << "  " << imag[i] << std::endl ;

            std::cout << std::endl << std::endl;
            delete[] real;
            delete[] imag;
            delete[] real1;
            delete[] imag1;
        }
    }

    std::cout << "Test24 end " << std::endl << std::endl;

}

/*
  We try to evaluate the polynomium:  1 + X + 3X^3  mod 11
  for a geometric  sequence  e, er, er^2, er^3, er^4...
  lets begin with e = 1 and r = 2
                             1,2,4,8,16

   some basic numbers  (all is modulo 11)

     a      a^2        a^-1     a^-3    a^-6
     0      0           n.a     n.a.    n.a
     1      1           1       1       1
     2      4           6       7       5
     3      9           4       9       4
     4      5           3       5       3
     5      3           9       3       9
     6      3           2       8       9
     7      5           8       4       5
     8      9           7       2       4
     9      4           5       4       5
     10     1           10      10      1



   y0 = a0 * r^-0    = 1
   y1 = a1 * r^-1    = 6
   y2 = a2 * r^-3    = 0
   y3 = a3 * r^-6    = 5
*/
/*
void test25()
{
    ResidueClass RC;

    PrimeFactorDFT fft;

    factorSeq  factors;

    fft.CalcFactors(1000, factors, 4);

    std::string  text = " 83264826X^10 - 23917X^2 + 27713X  + 9212937";

    RC.Setup((char *)   "100000000", 100000);

    PolynomialasResidue  PR(text, RC);

    std::cout << "Test25 end " << std::endl << std::endl;

}
*/

int main()
{
    //test1();
    //test2DNA();
    //test3Convolution();
    //test4();
    //test5perf();
    //test6Calc();
    //test7Calc();
    //test8Calc();
    //test9Calc();
    //test10Calc();
    //test11CalcGCD();
    //test12CalcSmall();
    //test13QuotientReminder();
    //test14Exp();
    //test15Jacobi();
    //test17();
    //test16MillerRabin((char*)"2147483647");// Mersenne Prime
    //test16MillerRabin((char *) "1228467");
    //test16MillerRabin((char*)"333228469");
    //test16MillerRabin((char*)"19777122847");
    //test18();
    //Factoring((char*)"2147483649");
    //test19();
    test20();
    //test21();
    //test22();
    //test23();
    //test24();
    //test25();
    std::cout << "Done !\n";
}
