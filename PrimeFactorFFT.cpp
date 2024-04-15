/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/
// for performance measurement

#include "PrimeFactorDFT.h"

#include "PrimeFactorDFT.h"
#include <iostream>
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

    std::cout << "Test4 begin " << std::endl;

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
    std::cout << "Test4 end " << std::endl << std::endl;

}


int main()
{
    test1();
    test2DNA();
    test3Convolution();
    test4();
    std::cout << "Done !\n";
}
