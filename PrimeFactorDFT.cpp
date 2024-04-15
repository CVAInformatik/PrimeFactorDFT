/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.

Inspiration originally from

                Clive Temperton
                Journal of Computational Physics 58, 283-299 (1985)
                "Implementation of a Self-Sorting In-Place Prime Factor FFT Algorithm."

         And the work by

                    C. Sidney Burrus
                    Ivan Selesnick   at RICE University

*/


#include <iostream>
#include <list>
#include "PrimeFactorDFT.h"

s64 PrimeFactorDFT::ValidateFactors(factorSeq& _factors)
{
    int f2 = 0;
    int f3 = 0;
	int f5 = 0;
	int f7 = 0;
	int f11 = 0;
	int f13 = 0;
	int f17 = 0;
	int f19 = 0;
	int f31 = 0;
	int funknown = 0;

	for (factorSeq::const_iterator it = _factors.begin();
		it != _factors.end();
		it++)
		switch (*it)
		{
        case 2: f2++; break;
        case 3: f3++; break;
		case 5: f5++; break;
		case 7: f7++; break;
		case 11: f11++; break;
		case 13: f13++; break;
		case 17: f17++; break;
		case 19: f19++; break;
		case 31: f31++; break;
		default: funknown++; break;
		};

	if (0 == (f2 + f3 + f5 + f7 + f11 + f13 + f17 + f19 + f31))
		return 0;
	if (funknown > 0)
		return -1;
	if ((f2 > 1) || (f3 > 1) || (f5 > 1) || (f7 > 1) || (f11 > 1) || (f13 > 1)
		|| (f17 > 1) || (f19 > 1) || (f31 > 1))
		return -2;

	s64 length = 1;
    if (f2 == 1) length *= 2;
    if (f3 == 1) length *= 3;
	if (f5 == 1) length *= 5;
	if (f7 == 1) length *= 7;
	if (f11 == 1) length *= 11;
	if (f13 == 1) length *= 13;
	if (f17 == 1) length *= 17;
	if (f19 == 1) length *= 19;
	if (f31 == 1) length *= 31;
	return length;
}

void PrimeFactorDFT::InitDFT(factorSeq& _factors, std::vector<BasicDFT*> &_DFTs)
{
	for (std::size_t i = 0;i < _factors.size();i++)
	{
		std::vector<s64>  indices;
		InitIndices(indices, _factors[i], (int) state);
		BasicDFT* t;
		switch (_factors[i])
		{
        case 2:  t = (BasicDFT*) new   DFT2(Rotations[i], state / 2, indices); 	_DFTs.push_back(t); break;
        case 3:  t = (BasicDFT*) new   DFT3(Rotations[i], state / 3, indices); 	_DFTs.push_back(t); break;
		case 5:  t = (BasicDFT*) new   DFT5(Rotations[i], state / 5, indices); 	_DFTs.push_back(t); break;
		case 7:  t = (BasicDFT*) new   DFT7(Rotations[i], state / 7, indices); 	_DFTs.push_back(t); break;
		case 11: t = (BasicDFT*) new  DFT11(Rotations[i], state / 11, indices); 	_DFTs.push_back(t); break;
		case 13: t = (BasicDFT*) new  DFT13(Rotations[i], state / 13, indices); 	_DFTs.push_back(t); break;
		case 17: t = (BasicDFT*) new  DFT17(Rotations[i], state / 17, indices); 	_DFTs.push_back(t); break;
		case 19: t = (BasicDFT*) new  DFT19(Rotations[i], state / 19, indices); 	_DFTs.push_back(t); break;
		case 31: t = (BasicDFT*) new  DFT31(Rotations[i], state / 31, indices); 	_DFTs.push_back(t); break;
		default: std::cout << "PFADFT::PFADT something is wrong here, Factorlist[" << i << "]= " << state << std::endl;
		}
	}

}

void PrimeFactorDFT::CleanUpDFT( std::vector<BasicDFT*> &_DTFs)
{}

void PrimeFactorDFT::forwardFFT(Data* real, Data *imag)
{
	for (std::vector<BasicDFT*>::const_iterator it = DFTs.begin();it != DFTs.end();it++)
	{
		(*it)->Evaluate(real, imag);
	}
};
void PrimeFactorDFT::InverseFFT(Data* real, Data *imag)
{
	for (std::vector<BasicDFT*>::const_iterator it = DFTs.begin();it != DFTs.end();it++)
	{
		(*it)->Evaluate(imag, real);
	}

};
void PrimeFactorDFT::ScaledInverseFFT(Data* real, Data *imag)
{
	for (std::vector<BasicDFT*>::const_iterator it = DFTs.begin();it != DFTs.end();it++)
	{
		(*it)->Evaluate(imag, real);
	}
	for (uint i = 0; i < state; i++)
	{
		real[i] /= state;
		imag[i] /= state;
	}
};


int PrimeFactorDFT::CalcFactors(uint length, factorSeq& _factors, int factorCount)
{
    std::list<unsigned int> lengthList;


    for(int ix = 0; ix < 512; ix++) {

        int hw = 0;
        if (factorCount) {
            int it = ix;
            while (it) {
                it = it & (it - 1);
                hw++;
            }
        }
        if (factorCount && ( hw > factorCount)) continue;

        uint tlength = 1;
        if (ix & 1) tlength *= 2;
        if (ix & 2) tlength *= 3;
        if (ix & 4) tlength *= 5;
        if (ix & 8) tlength *= 7;
        if (ix & 16) tlength *= 11;
        if (ix & 32) tlength *= 13;
        if (ix & 64) tlength *= 17;
        if (ix & 128) tlength *= 19;
        if (ix & 256) tlength *= 31;
        lengthList.push_back(tlength);
    }
    lengthList.sort();

    uint actualLength = 0;
    _factors.clear();

    for (std::list<unsigned int>::const_iterator ibegin = lengthList.begin(); 1 ; ibegin++)
        if (ibegin == lengthList.end())  return -1;        
        else if (*ibegin >= length)
        {
            actualLength = *ibegin;
            break;
        }

    if ((actualLength % 2 )== 0) _factors.push_back(2);
    if ((actualLength % 3) == 0) _factors.push_back(3);
    if ((actualLength % 5) == 0) _factors.push_back(5);
    if ((actualLength % 7) == 0) _factors.push_back(7);
    if ((actualLength % 11) == 0) _factors.push_back(11);
    if ((actualLength % 13) == 0) _factors.push_back(13);
    if ((actualLength % 17) == 0) _factors.push_back(17);
    if ((actualLength % 19) == 0) _factors.push_back(19);
    if ((actualLength % 31) == 0) _factors.push_back(31);

    return actualLength;
}


/*
	Code from Temperton

	Algebraically the rotation  for a factor is the multiplicative inverse  of   N/factor  mod factor
	for N = 5040 	  9 x 560   560  mod 9  ==  2  and  2 x 5  = 1 mod 9 so rotation is 5  for the length 9 DFT
					 16 x 315   315  mod 16 == 11  and  3 x 11 = 1 mod 16 so rotation is 3  for the length 16 DFT
					  5 x 1008  1008 mod 5  ==  3  and  2 x 3  = 1 mod 5  so rotation is 2  for the length 5 DFT
					  7 x 720   720  mod 7  ==  6  and  6 x 6  = 1 mod 7  so rotation is 6  for the length 7 DFT

	The code below is (almost) the original Fortran IV code from Temperton. (Fortran indexes from 1 !)

	*/
void PrimeFactorDFT::InitRotations()
{
	int N = 1;

	Rotations.clear();
	for (factorSeq::const_iterator cit = factors.begin(); cit != factors.end(); cit++)  N = *cit * N;

	for (factorSeq::const_iterator cit = factors.begin(); cit != factors.end(); cit++) {

		int IFAC = (int)*cit;
		int M = N / IFAC;
		int MU = 0;
		int MM = 0;
		for (int J = 1; J < IFAC; J++)
		{
			MU = J;
			MM = J * M;
			if (MM % IFAC == 1)
			{
				Rotations.push_back(MU);
				break;
			}
		}
	}

}

void PrimeFactorDFT::InitIndices(std::vector<s64>& indices, int fftlength, s64 length)
{
    indices.clear();
    indices.resize(fftlength);

    s64 offset = 0;

    while (offset < length)
    {
        indices[offset % fftlength] = offset;
        offset += (length / fftlength);
    }
}

void DFT2::Evaluate(Data *real, Data *imag)
{
    std::vector<s64> ind = indices;

    for (int i = 0; i < count; i++)
    {
        Data t1real = real[ind[0]] + real[ind[1]];
        Data t1imag = imag[ind[0]] + imag[ind[1]];
        Data t2real = real[ind[0]] - real[ind[1]];
        Data t2imag = imag[ind[0]] - imag[ind[1]];
        real[ind[0]] = t1real;
        imag[ind[0]] = t1imag;
        real[ind[1]] = t2real;
        imag[ind[1]] = t2imag;
        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}

#undef FFTLENGTH
#define FFTLENGTH 3

void DFT3::Evaluate(Data* real, Data *imag)
{
    Data real_x[FFTLENGTH];
    Data imag_x[FFTLENGTH];



    Data real_y[FFTLENGTH];
    Data imag_y[FFTLENGTH];

    Data real_v[FFTLENGTH];
    Data imag_v[FFTLENGTH];

    std::vector<s64> ind = indices;

    for (int i = 0; i < count; i++)
    {
        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real[ind[px]];
            imag_x[px] = imag[ind[px]];
        }
        /* input permuation is void */

        //
        // DFT length 3
        //
         /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[2] = real_x[1] - real_x[2];
        /* RED */       imag_v[2] = imag_x[1] - imag_x[2];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULim  */       real_y[2] = -1 * imag_x[2] * u[1];
        /* MULim  */       imag_y[2] = real_x[2] * u[1];
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* KRED */       /* tKRED entry */
          /* tRED */       /*  tRED entry p 2 a 1 c 1  */
        /* tRED */       real_v[2] = real_y[1];
        /* tRED */       imag_v[2] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[2];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[2];
        /* tRED */       real_v[2] -= real_y[2];
        /* tRED */       imag_v[2] -= imag_y[2];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       /*  tRED exit  */
      /* KRED */       /*  tKRED exit */

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[active_op[px]];
            imag_x[px] = imag_y[active_op[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }


        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}


#undef FFTLENGTH
#define FFTLENGTH 5

void DFT5::Evaluate(Data* real, Data* imag)
{
    Data real_x[FFTLENGTH];
    Data imag_x[FFTLENGTH];

    Data real_v[FFTLENGTH];
    Data imag_v[FFTLENGTH];
    Data real_y[FFTLENGTH];
    Data imag_y[FFTLENGTH];
    Data real_t, imag_t;

    std::vector<s64> ind = indices;


    for (int i = 0; i < count; i++)
    {
        for (int px = 0; px < FFTLENGTH; px++) {
            real_y[px] = real[ind[px]];
            imag_y[px] = imag[ind[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[ip[px]];
            imag_x[px] = imag_y[ip[px]];
        }

        //
        // DFT length 5
        //
  /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 2  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[3];
        /* RED */       imag_v[1] += imag_x[3];
        /* RED */       real_v[3] = real_x[1] - real_x[3];
        /* RED */       imag_v[3] = imag_x[1] - imag_x[3];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[4];
        /* RED */       imag_v[2] += imag_x[4];
        /* RED */       real_v[4] = real_x[2] - real_x[4];
        /* RED */       imag_v[4] = imag_x[2] - imag_x[4];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 2  a 1 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[2] = real_x[1] - real_x[2];
        /* RED */       imag_v[2] = imag_x[1] - imag_x[2];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULre  */       real_y[2] = real_x[2] * u[1];
        /* MULre  */       imag_y[2] = imag_x[2] * u[1];
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[3];
        /* D2  */       imag_v[0] = imag_x[3];
        /* D2  */       real_v[1] = real_x[4];
        /* D2  */       imag_v[1] = imag_x[4];
        /* D2  */       real_v[2] = real_x[3] + real_x[4];
        /* D2  */       imag_v[2] = imag_x[3] + imag_x[4];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[2];
        /* IMAG  */       imag_t = real_v[0] * u[2];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[3];
        /* IMAG  */       imag_t = real_v[1] * u[3];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[4];
        /* IMAG  */       imag_t = real_v[2] * u[4];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[3] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[3] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[4] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[4] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* tKRED */       /* tKRED entry */
          /* tRED */       /*  tRED entry p 2 a 1 c 1  */
        /* tRED */       real_v[2] = real_y[1];
        /* tRED */       imag_v[2] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[2];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[2];
        /* tRED */       real_v[2] -= real_y[2];
        /* tRED */       imag_v[2] -= imag_y[2];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 2  */
        /* tRED */       real_v[3] = real_y[1];
        /* tRED */       imag_v[3] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[3];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[3];
        /* tRED */       real_v[3] -= real_y[3];
        /* tRED */       imag_v[3] -= imag_y[3];
        /* tRED */       real_v[4] = real_y[2];
        /* tRED */       imag_v[4] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[4];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[4];
        /* tRED */       real_v[4] -= real_y[4];
        /* tRED */       imag_v[4] -= imag_y[4];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       /*  tRED exit  */
      /* tKRED */       /*  tKRED exit */


        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[active_op[px]];
            imag_x[px] = imag_y[active_op[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }


        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}

#undef FFTLENGTH
#define FFTLENGTH 7

void DFT7::Evaluate(Data* real, Data* imag)
{
    Data real_x[FFTLENGTH];
    Data imag_x[FFTLENGTH];

    Data real_v[FFTLENGTH];
    Data imag_v[FFTLENGTH];
    Data real_v1[FFTLENGTH];
    Data imag_v1[FFTLENGTH];
    Data real_y[FFTLENGTH];
    Data imag_y[FFTLENGTH];
    Data real_t, imag_t;

    std::vector<s64> ind = indices;


    for (int i = 0; i < count; i++)
    {
        for (int px = 0; px < FFTLENGTH; px++) {
            real_y[px] = real[ind[px]];
            imag_y[px] = imag[ind[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[ip[px]];
            imag_x[px] = imag_y[ip[px]];
        }

        //
        // DFT length 7
        //

  /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 3  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[4];
        /* RED */       imag_v[1] += imag_x[4];
        /* RED */       real_v[4] = real_x[1] - real_x[4];
        /* RED */       imag_v[4] = imag_x[1] - imag_x[4];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[5];
        /* RED */       imag_v[2] += imag_x[5];
        /* RED */       real_v[5] = real_x[2] - real_x[5];
        /* RED */       imag_v[5] = imag_x[2] - imag_x[5];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[6];
        /* RED */       imag_v[3] += imag_x[6];
        /* RED */       real_v[6] = real_x[3] - real_x[6];
        /* RED */       imag_v[6] = imag_x[3] - imag_x[6];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 3  a 2 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[3] = real_x[1] - real_x[3];
        /* RED */       imag_v[3] = imag_x[1] - imag_x[3];
        /* RED */       real_v[1] += real_x[3];
        /* RED */       imag_v[1] += imag_x[3];
        /* RED */       real_v[4] = real_x[2] - real_x[3];
        /* RED */       imag_v[4] = imag_x[2] - imag_x[3];
        /* RED */       real_v[2] = real_x[4];
        /* RED */       imag_v[2] = imag_x[4];
        /* RED */       real_v[2] += real_x[5];
        /* RED */       imag_v[2] += imag_x[5];
        /* RED */       real_v[5] = real_x[4] - real_x[6];
        /* RED */       imag_v[5] = imag_x[4] - imag_x[6];
        /* RED */       real_v[2] += real_x[6];
        /* RED */       imag_v[2] += imag_x[6];
        /* RED */       real_v[6] = real_x[5] - real_x[6];
        /* RED */       imag_v[6] = imag_x[5] - imag_x[6];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULim  */       real_y[2] = -1 * imag_x[2] * u[1];
        /* MULim  */       imag_y[2] = real_x[2] * u[1];
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[3];
        /* D2  */       imag_v[0] = imag_x[3];
        /* D2  */       real_v[1] = real_x[4];
        /* D2  */       imag_v[1] = imag_x[4];
        /* D2  */       real_v[2] = real_x[3] + real_x[4];
        /* D2  */       imag_v[2] = imag_x[3] + imag_x[4];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[2];
        /* REAL  */       imag_v[0] = imag_v[0] * u[2];
        /* REAL  */       real_v[1] = real_v[1] * u[3];
        /* REAL  */       imag_v[1] = imag_v[1] * u[3];
        /* REAL  */       real_v[2] = real_v[2] * u[4];
        /* REAL  */       imag_v[2] = imag_v[2] * u[4];
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[3] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[3] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[4] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[4] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v1[0] = real_v[5];
        /* D2  */       imag_v1[0] = imag_v[5];
        /* D2  */       real_v1[1] = real_v[6];
        /* D2  */       imag_v1[1] = imag_v[6];
        /* D2  */       real_v1[2] = real_v[5] + real_v[6];
        /* D2  */       imag_v1[2] = imag_v[5] + imag_v[6];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v1[0] * u[5];
        /* IMAG  */       imag_t = real_v1[0] * u[5];
        /* IMAG  */       real_v1[0] = real_t;
        /* IMAG  */       imag_v1[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[1] * u[6];
        /* IMAG  */       imag_t = real_v1[1] * u[6];
        /* IMAG  */       real_v1[1] = real_t;
        /* IMAG  */       imag_v1[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[2] * u[7];
        /* IMAG  */       imag_t = real_v1[2] * u[7];
        /* IMAG  */       real_v1[2] = real_t;
        /* IMAG  */       imag_v1[2] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[5] = real_v1[0] + real_v1[2];
        /* D2t  */       imag_y[5] = imag_v1[0] + imag_v1[2];
        /* D2t  */       real_y[6] = real_v1[1] + real_v1[2];
        /* D2t  */       imag_y[6] = imag_v1[1] + imag_v1[2];
        /* ID2It */     /* Exit */
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* tKRED */       /* tKRED entry */
          /* tRED */       /*  tRED entry p 3 a 2 c 1  */
        /* tRED */       real_v[3] = real_y[1];
        /* tRED */       imag_v[3] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[3];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[3];
        /* tRED */       real_v[3] -= real_y[3];
        /* tRED */       imag_v[3] -= imag_y[3];
        /* tRED */       real_v[2] = real_y[1] + real_y[4];
        /* tRED */       imag_v[2] = imag_y[1] + imag_y[4];
        /* tRED */       real_v[3] -= real_y[4];
        /* tRED */       imag_v[3] -= imag_y[4];
        /* tRED */       real_v[6] = real_y[2];
        /* tRED */       imag_v[6] = imag_y[2];
        /* tRED */       real_v[4] = real_y[2] + real_y[5];
        /* tRED */       imag_v[4] = imag_y[2] + imag_y[5];
        /* tRED */       real_v[6] -= real_y[5];
        /* tRED */       imag_v[6] -= imag_y[5];
        /* tRED */       real_v[5] = real_y[2] + real_y[6];
        /* tRED */       imag_v[5] = imag_y[2] + imag_y[6];
        /* tRED */       real_v[6] -= real_y[6];
        /* tRED */       imag_v[6] -= imag_y[6];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 3  */
        /* tRED */       real_v[4] = real_y[1];
        /* tRED */       imag_v[4] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[4];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[4];
        /* tRED */       real_v[4] -= real_y[4];
        /* tRED */       imag_v[4] -= imag_y[4];
        /* tRED */       real_v[5] = real_y[2];
        /* tRED */       imag_v[5] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[5];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[5];
        /* tRED */       real_v[5] -= real_y[5];
        /* tRED */       imag_v[5] -= imag_y[5];
        /* tRED */       real_v[6] = real_y[3];
        /* tRED */       imag_v[6] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[6];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[6];
        /* tRED */       real_v[6] -= real_y[6];
        /* tRED */       imag_v[6] -= imag_y[6];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       /*  tRED exit  */
      /* tKRED */       /*  tKRED exit */

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[active_op[px]];
            imag_x[px] = imag_y[active_op[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }


        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}



#undef FFTLENGTH
#define FFTLENGTH 11

void DFT11::Evaluate(Data* real, Data* imag)
{
    std::vector<s64> ind = indices;


	Data real_x[FFTLENGTH];
	Data imag_x[FFTLENGTH];

	Data real_v[FFTLENGTH];
	Data imag_v[FFTLENGTH];
	Data real_v1[FFTLENGTH];
	Data imag_v1[FFTLENGTH];
	Data real_y[FFTLENGTH];
	Data imag_y[FFTLENGTH];
    Data real_t, imag_t;


	for (int i = 0; i < count; i++)
	{

		/// OBS OBS mangler rotatation lige nu 

        for (int px = 0; px < FFTLENGTH; px++) {
            real_y[px] = real[ind[px]];
            imag_y[px] = imag[ind[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[ip[px]];
            imag_x[px] = imag_y[ip[px]];
        }

		//
		// DFT length 11
		//
		//EVALUATEDFT(11);
  /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 5  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[6];
        /* RED */       imag_v[1] += imag_x[6];
        /* RED */       real_v[6] = real_x[1] - real_x[6];
        /* RED */       imag_v[6] = imag_x[1] - imag_x[6];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[7];
        /* RED */       imag_v[2] += imag_x[7];
        /* RED */       real_v[7] = real_x[2] - real_x[7];
        /* RED */       imag_v[7] = imag_x[2] - imag_x[7];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[8];
        /* RED */       imag_v[3] += imag_x[8];
        /* RED */       real_v[8] = real_x[3] - real_x[8];
        /* RED */       imag_v[8] = imag_x[3] - imag_x[8];
        /* RED */       real_v[4] = real_x[4];
        /* RED */       imag_v[4] = imag_x[4];
        /* RED */       real_v[4] += real_x[9];
        /* RED */       imag_v[4] += imag_x[9];
        /* RED */       real_v[9] = real_x[4] - real_x[9];
        /* RED */       imag_v[9] = imag_x[4] - imag_x[9];
        /* RED */       real_v[5] = real_x[5];
        /* RED */       imag_v[5] = imag_x[5];
        /* RED */       real_v[5] += real_x[10];
        /* RED */       imag_v[5] += imag_x[10];
        /* RED */       real_v[10] = real_x[5] - real_x[10];
        /* RED */       imag_v[10] = imag_x[5] - imag_x[10];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 5  a 2 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[3] = real_x[1] - real_x[5];
        /* RED */       imag_v[3] = imag_x[1] - imag_x[5];
        /* RED */       real_v[1] += real_x[3];
        /* RED */       imag_v[1] += imag_x[3];
        /* RED */       real_v[4] = real_x[2] - real_x[5];
        /* RED */       imag_v[4] = imag_x[2] - imag_x[5];
        /* RED */       real_v[1] += real_x[4];
        /* RED */       imag_v[1] += imag_x[4];
        /* RED */       real_v[5] = real_x[3] - real_x[5];
        /* RED */       imag_v[5] = imag_x[3] - imag_x[5];
        /* RED */       real_v[1] += real_x[5];
        /* RED */       imag_v[1] += imag_x[5];
        /* RED */       real_v[6] = real_x[4] - real_x[5];
        /* RED */       imag_v[6] = imag_x[4] - imag_x[5];
        /* RED */       real_v[2] = real_x[6];
        /* RED */       imag_v[2] = imag_x[6];
        /* RED */       real_v[2] += real_x[7];
        /* RED */       imag_v[2] += imag_x[7];
        /* RED */       real_v[7] = real_x[6] - real_x[10];
        /* RED */       imag_v[7] = imag_x[6] - imag_x[10];
        /* RED */       real_v[2] += real_x[8];
        /* RED */       imag_v[2] += imag_x[8];
        /* RED */       real_v[8] = real_x[7] - real_x[10];
        /* RED */       imag_v[8] = imag_x[7] - imag_x[10];
        /* RED */       real_v[2] += real_x[9];
        /* RED */       imag_v[2] += imag_x[9];
        /* RED */       real_v[9] = real_x[8] - real_x[10];
        /* RED */       imag_v[9] = imag_x[8] - imag_x[10];
        /* RED */       real_v[2] += real_x[10];
        /* RED */       imag_v[2] += imag_x[10];
        /* RED */       real_v[10] = real_x[9] - real_x[10];
        /* RED */       imag_v[10] = imag_x[9] - imag_x[10];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULim  */       real_y[2] = -1 * imag_x[2] * u[1];
        /* MULim  */       imag_y[2] = real_x[2] * u[1];
        /* ID2I */     /* Entry m = 1  n = 2 */
        /* D2  */       real_v[0] = real_x[3];
        /* D2  */       imag_v[0] = imag_x[3];
        /* D2  */       real_v[2] = real_x[5];
        /* D2  */       imag_v[2] = imag_x[5];
        /* D2  */       real_v[4] = real_x[3] + real_x[5];
        /* D2  */       imag_v[4] = imag_x[3] + imag_x[5];
        /* D2  */       real_v[1] = real_x[4];
        /* D2  */       imag_v[1] = imag_x[4];
        /* D2  */       real_v[3] = real_x[6];
        /* D2  */       imag_v[3] = imag_x[6];
        /* D2  */       real_v[5] = real_x[4] + real_x[6];
        /* D2  */       imag_v[5] = imag_x[4] + imag_x[6];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 1 */
        /* D2  */       real_v1[0] = real_v[0];
        /* D2  */       imag_v1[0] = imag_v[0];
        /* D2  */       real_v1[1] = real_v[1];
        /* D2  */       imag_v1[1] = imag_v[1];
        /* D2  */       real_v1[2] = real_v[0] + real_v[1];
        /* D2  */       imag_v1[2] = imag_v[0] + imag_v[1];
        /* D2  */       real_v1[3] = real_v[2];
        /* D2  */       imag_v1[3] = imag_v[2];
        /* D2  */       real_v1[4] = real_v[3];
        /* D2  */       imag_v1[4] = imag_v[3];
        /* D2  */       real_v1[5] = real_v[2] + real_v[3];
        /* D2  */       imag_v1[5] = imag_v[2] + imag_v[3];
        /* D2  */       real_v1[6] = real_v[4];
        /* D2  */       imag_v1[6] = imag_v[4];
        /* D2  */       real_v1[7] = real_v[5];
        /* D2  */       imag_v1[7] = imag_v[5];
        /* D2  */       real_v1[8] = real_v[4] + real_v[5];
        /* D2  */       imag_v1[8] = imag_v[4] + imag_v[5];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v1[0] = real_v1[0] * u[2];
        /* REAL  */       imag_v1[0] = imag_v1[0] * u[2];
        /* REAL  */       real_v1[1] = real_v1[1] * u[3];
        /* REAL  */       imag_v1[1] = imag_v1[1] * u[3];
        /* REAL  */       real_v1[2] = real_v1[2] * u[4];
        /* REAL  */       imag_v1[2] = imag_v1[2] * u[4];
        /* REAL  */       real_v1[3] = real_v1[3] * u[5];
        /* REAL  */       imag_v1[3] = imag_v1[3] * u[5];
        /* REAL  */       real_v1[4] = real_v1[4] * u[6];
        /* REAL  */       imag_v1[4] = imag_v1[4] * u[6];
        /* REAL  */       real_v1[5] = real_v1[5] * u[7];
        /* REAL  */       imag_v1[5] = imag_v1[5] * u[7];
        /* REAL  */       real_v1[6] = real_v1[6] * u[8];
        /* REAL  */       imag_v1[6] = imag_v1[6] * u[8];
        /* REAL  */       real_v1[7] = real_v1[7] * u[9];
        /* REAL  */       imag_v1[7] = imag_v1[7] * u[9];
        /* REAL  */       real_v1[8] = real_v1[8] * u[10];
        /* REAL  */       imag_v1[8] = imag_v1[8] * u[10];
        /* ID2It */     /* Entry m = 1  n = 3  */
        /* D2t  */       real_v[0] = real_v1[0] + real_v1[6];
        /* D2t  */       imag_v[0] = imag_v1[0] + imag_v1[6];
        /* D2t  */       real_v[3] = real_v1[3] + real_v1[6];
        /* D2t  */       imag_v[3] = imag_v1[3] + imag_v1[6];
        /* D2t  */       real_v[1] = real_v1[1] + real_v1[7];
        /* D2t  */       imag_v[1] = imag_v1[1] + imag_v1[7];
        /* D2t  */       real_v[4] = real_v1[4] + real_v1[7];
        /* D2t  */       imag_v[4] = imag_v1[4] + imag_v1[7];
        /* D2t  */       real_v[2] = real_v1[2] + real_v1[8];
        /* D2t  */       imag_v[2] = imag_v1[2] + imag_v1[8];
        /* D2t  */       real_v[5] = real_v1[5] + real_v1[8];
        /* D2t  */       imag_v[5] = imag_v1[5] + imag_v1[8];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 1  */
        /* D2t  */       real_y[3] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[3] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[4] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[4] = imag_v[1] + imag_v[2];
        /* D2t  */       real_y[5] = real_v[3] + real_v[5];
        /* D2t  */       imag_y[5] = imag_v[3] + imag_v[5];
        /* D2t  */       real_y[6] = real_v[4] + real_v[5];
        /* D2t  */       imag_y[6] = imag_v[4] + imag_v[5];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 2 */
        /* D2  */       real_v[0] = real_x[7];
        /* D2  */       imag_v[0] = imag_x[7];
        /* D2  */       real_v[2] = real_x[9];
        /* D2  */       imag_v[2] = imag_x[9];
        /* D2  */       real_v[4] = real_x[7] + real_x[9];
        /* D2  */       imag_v[4] = imag_x[7] + imag_x[9];
        /* D2  */       real_v[1] = real_x[8];
        /* D2  */       imag_v[1] = imag_x[8];
        /* D2  */       real_v[3] = real_x[10];
        /* D2  */       imag_v[3] = imag_x[10];
        /* D2  */       real_v[5] = real_x[8] + real_x[10];
        /* D2  */       imag_v[5] = imag_x[8] + imag_x[10];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 1 */
        /* D2  */       real_v1[0] = real_v[0];
        /* D2  */       imag_v1[0] = imag_v[0];
        /* D2  */       real_v1[1] = real_v[1];
        /* D2  */       imag_v1[1] = imag_v[1];
        /* D2  */       real_v1[2] = real_v[0] + real_v[1];
        /* D2  */       imag_v1[2] = imag_v[0] + imag_v[1];
        /* D2  */       real_v1[3] = real_v[2];
        /* D2  */       imag_v1[3] = imag_v[2];
        /* D2  */       real_v1[4] = real_v[3];
        /* D2  */       imag_v1[4] = imag_v[3];
        /* D2  */       real_v1[5] = real_v[2] + real_v[3];
        /* D2  */       imag_v1[5] = imag_v[2] + imag_v[3];
        /* D2  */       real_v1[6] = real_v[4];
        /* D2  */       imag_v1[6] = imag_v[4];
        /* D2  */       real_v1[7] = real_v[5];
        /* D2  */       imag_v1[7] = imag_v[5];
        /* D2  */       real_v1[8] = real_v[4] + real_v[5];
        /* D2  */       imag_v1[8] = imag_v[4] + imag_v[5];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v1[0] * u[11];
        /* IMAG  */       imag_t = real_v1[0] * u[11];
        /* IMAG  */       real_v1[0] = real_t;
        /* IMAG  */       imag_v1[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[1] * u[12];
        /* IMAG  */       imag_t = real_v1[1] * u[12];
        /* IMAG  */       real_v1[1] = real_t;
        /* IMAG  */       imag_v1[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[2] * u[13];
        /* IMAG  */       imag_t = real_v1[2] * u[13];
        /* IMAG  */       real_v1[2] = real_t;
        /* IMAG  */       imag_v1[2] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[3] * u[14];
        /* IMAG  */       imag_t = real_v1[3] * u[14];
        /* IMAG  */       real_v1[3] = real_t;
        /* IMAG  */       imag_v1[3] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[4] * u[15];
        /* IMAG  */       imag_t = real_v1[4] * u[15];
        /* IMAG  */       real_v1[4] = real_t;
        /* IMAG  */       imag_v1[4] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[5] * u[16];
        /* IMAG  */       imag_t = real_v1[5] * u[16];
        /* IMAG  */       real_v1[5] = real_t;
        /* IMAG  */       imag_v1[5] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[6] * u[17];
        /* IMAG  */       imag_t = real_v1[6] * u[17];
        /* IMAG  */       real_v1[6] = real_t;
        /* IMAG  */       imag_v1[6] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[7] * u[18];
        /* IMAG  */       imag_t = real_v1[7] * u[18];
        /* IMAG  */       real_v1[7] = real_t;
        /* IMAG  */       imag_v1[7] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[8] * u[19];
        /* IMAG  */       imag_t = real_v1[8] * u[19];
        /* IMAG  */       real_v1[8] = real_t;
        /* IMAG  */       imag_v1[8] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 3  */
        /* D2t  */       real_v[0] = real_v1[0] + real_v1[6];
        /* D2t  */       imag_v[0] = imag_v1[0] + imag_v1[6];
        /* D2t  */       real_v[3] = real_v1[3] + real_v1[6];
        /* D2t  */       imag_v[3] = imag_v1[3] + imag_v1[6];
        /* D2t  */       real_v[1] = real_v1[1] + real_v1[7];
        /* D2t  */       imag_v[1] = imag_v1[1] + imag_v1[7];
        /* D2t  */       real_v[4] = real_v1[4] + real_v1[7];
        /* D2t  */       imag_v[4] = imag_v1[4] + imag_v1[7];
        /* D2t  */       real_v[2] = real_v1[2] + real_v1[8];
        /* D2t  */       imag_v[2] = imag_v1[2] + imag_v1[8];
        /* D2t  */       real_v[5] = real_v1[5] + real_v1[8];
        /* D2t  */       imag_v[5] = imag_v1[5] + imag_v1[8];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 1  */
        /* D2t  */       real_y[7] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[7] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[8] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[8] = imag_v[1] + imag_v[2];
        /* D2t  */       real_y[9] = real_v[3] + real_v[5];
        /* D2t  */       imag_y[9] = imag_v[3] + imag_v[5];
        /* D2t  */       real_y[10] = real_v[4] + real_v[5];
        /* D2t  */       imag_y[10] = imag_v[4] + imag_v[5];
        /* ID2It */     /* Exit */
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* tKRED */       /* tKRED entry */
          /* tRED */       /*  tRED entry p 5 a 2 c 1  */
        /* tRED */       real_v[5] = real_y[1];
        /* tRED */       imag_v[5] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[3];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[3];
        /* tRED */       real_v[5] -= real_y[3];
        /* tRED */       imag_v[5] -= imag_y[3];
        /* tRED */       real_v[2] = real_y[1] + real_y[4];
        /* tRED */       imag_v[2] = imag_y[1] + imag_y[4];
        /* tRED */       real_v[5] -= real_y[4];
        /* tRED */       imag_v[5] -= imag_y[4];
        /* tRED */       real_v[3] = real_y[1] + real_y[5];
        /* tRED */       imag_v[3] = imag_y[1] + imag_y[5];
        /* tRED */       real_v[5] -= real_y[5];
        /* tRED */       imag_v[5] -= imag_y[5];
        /* tRED */       real_v[4] = real_y[1] + real_y[6];
        /* tRED */       imag_v[4] = imag_y[1] + imag_y[6];
        /* tRED */       real_v[5] -= real_y[6];
        /* tRED */       imag_v[5] -= imag_y[6];
        /* tRED */       real_v[10] = real_y[2];
        /* tRED */       imag_v[10] = imag_y[2];
        /* tRED */       real_v[6] = real_y[2] + real_y[7];
        /* tRED */       imag_v[6] = imag_y[2] + imag_y[7];
        /* tRED */       real_v[10] -= real_y[7];
        /* tRED */       imag_v[10] -= imag_y[7];
        /* tRED */       real_v[7] = real_y[2] + real_y[8];
        /* tRED */       imag_v[7] = imag_y[2] + imag_y[8];
        /* tRED */       real_v[10] -= real_y[8];
        /* tRED */       imag_v[10] -= imag_y[8];
        /* tRED */       real_v[8] = real_y[2] + real_y[9];
        /* tRED */       imag_v[8] = imag_y[2] + imag_y[9];
        /* tRED */       real_v[10] -= real_y[9];
        /* tRED */       imag_v[10] -= imag_y[9];
        /* tRED */       real_v[9] = real_y[2] + real_y[10];
        /* tRED */       imag_v[9] = imag_y[2] + imag_y[10];
        /* tRED */       real_v[10] -= real_y[10];
        /* tRED */       imag_v[10] -= imag_y[10];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 5  */
        /* tRED */       real_v[6] = real_y[1];
        /* tRED */       imag_v[6] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[6];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[6];
        /* tRED */       real_v[6] -= real_y[6];
        /* tRED */       imag_v[6] -= imag_y[6];
        /* tRED */       real_v[7] = real_y[2];
        /* tRED */       imag_v[7] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[7];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[7];
        /* tRED */       real_v[7] -= real_y[7];
        /* tRED */       imag_v[7] -= imag_y[7];
        /* tRED */       real_v[8] = real_y[3];
        /* tRED */       imag_v[8] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[8];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[8];
        /* tRED */       real_v[8] -= real_y[8];
        /* tRED */       imag_v[8] -= imag_y[8];
        /* tRED */       real_v[9] = real_y[4];
        /* tRED */       imag_v[9] = imag_y[4];
        /* tRED */       real_v[4] = real_y[4] + real_y[9];
        /* tRED */       imag_v[4] = imag_y[4] + imag_y[9];
        /* tRED */       real_v[9] -= real_y[9];
        /* tRED */       imag_v[9] -= imag_y[9];
        /* tRED */       real_v[10] = real_y[5];
        /* tRED */       imag_v[10] = imag_y[5];
        /* tRED */       real_v[5] = real_y[5] + real_y[10];
        /* tRED */       imag_v[5] = imag_y[5] + imag_y[10];
        /* tRED */       real_v[10] -= real_y[10];
        /* tRED */       imag_v[10] -= imag_y[10];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       /*  tRED exit  */
      /* tKRED */       /*  tKRED exit */





		//data[op]		
        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[active_op[px]];
            imag_x[px] = imag_y[active_op[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }


		//
		//  CRT mapping.
		//
		IncIndices(ind);
	}
}

#undef FFTLENGTH
#define FFTLENGTH 13



void DFT13::Evaluate(Data* real, Data* imag)
{
    std::vector<s64> ind = indices;

    Data real_x[FFTLENGTH];
    Data imag_x[FFTLENGTH];
    Data real_v[FFTLENGTH];
    Data imag_v[FFTLENGTH];
    Data real_v1[18];
    Data imag_v1[18];
    Data real_y[FFTLENGTH];
    Data imag_y[FFTLENGTH];
    Data real_t, imag_t;

    for (int i = 0; i < count; i++)
    {

        for (int px = 0; px < FFTLENGTH; px++) {
            real_y[px] = real[ind[px]];
            imag_y[px] = imag[ind[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[ip[px]];
            imag_x[px] = imag_y[ip[px]];
        }

        //
        // DFT length 13
        //
  /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 6  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[7];
        /* RED */       imag_v[1] += imag_x[7];
        /* RED */       real_v[7] = real_x[1] - real_x[7];
        /* RED */       imag_v[7] = imag_x[1] - imag_x[7];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[8];
        /* RED */       imag_v[2] += imag_x[8];
        /* RED */       real_v[8] = real_x[2] - real_x[8];
        /* RED */       imag_v[8] = imag_x[2] - imag_x[8];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[9];
        /* RED */       imag_v[3] += imag_x[9];
        /* RED */       real_v[9] = real_x[3] - real_x[9];
        /* RED */       imag_v[9] = imag_x[3] - imag_x[9];
        /* RED */       real_v[4] = real_x[4];
        /* RED */       imag_v[4] = imag_x[4];
        /* RED */       real_v[4] += real_x[10];
        /* RED */       imag_v[4] += imag_x[10];
        /* RED */       real_v[10] = real_x[4] - real_x[10];
        /* RED */       imag_v[10] = imag_x[4] - imag_x[10];
        /* RED */       real_v[5] = real_x[5];
        /* RED */       imag_v[5] = imag_x[5];
        /* RED */       real_v[5] += real_x[11];
        /* RED */       imag_v[5] += imag_x[11];
        /* RED */       real_v[11] = real_x[5] - real_x[11];
        /* RED */       imag_v[11] = imag_x[5] - imag_x[11];
        /* RED */       real_v[6] = real_x[6];
        /* RED */       imag_v[6] = imag_x[6];
        /* RED */       real_v[6] += real_x[12];
        /* RED */       imag_v[6] += imag_x[12];
        /* RED */       real_v[12] = real_x[6] - real_x[12];
        /* RED */       imag_v[12] = imag_x[6] - imag_x[12];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 2  a 1 c 3  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[4];
        /* RED */       imag_v[1] += imag_x[4];
        /* RED */       real_v[4] = real_x[1] - real_x[4];
        /* RED */       imag_v[4] = imag_x[1] - imag_x[4];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[5];
        /* RED */       imag_v[2] += imag_x[5];
        /* RED */       real_v[5] = real_x[2] - real_x[5];
        /* RED */       imag_v[5] = imag_x[2] - imag_x[5];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[6];
        /* RED */       imag_v[3] += imag_x[6];
        /* RED */       real_v[6] = real_x[3] - real_x[6];
        /* RED */       imag_v[6] = imag_x[3] - imag_x[6];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 3  a 4 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[5] = real_x[1] - real_x[3];
        /* RED */       imag_v[5] = imag_x[1] - imag_x[3];
        /* RED */       real_v[1] += real_x[3];
        /* RED */       imag_v[1] += imag_x[3];
        /* RED */       real_v[6] = real_x[2] - real_x[3];
        /* RED */       imag_v[6] = imag_x[2] - imag_x[3];
        /* RED */       real_v[2] = real_x[4];
        /* RED */       imag_v[2] = imag_x[4];
        /* RED */       real_v[2] += real_x[5];
        /* RED */       imag_v[2] += imag_x[5];
        /* RED */       real_v[7] = real_x[4] - real_x[6];
        /* RED */       imag_v[7] = imag_x[4] - imag_x[6];
        /* RED */       real_v[2] += real_x[6];
        /* RED */       imag_v[2] += imag_x[6];
        /* RED */       real_v[8] = real_x[5] - real_x[6];
        /* RED */       imag_v[8] = imag_x[5] - imag_x[6];
        /* RED */       real_v[3] = real_x[7];
        /* RED */       imag_v[3] = imag_x[7];
        /* RED */       real_v[3] += real_x[8];
        /* RED */       imag_v[3] += imag_x[8];
        /* RED */       real_v[9] = real_x[7] - real_x[9];
        /* RED */       imag_v[9] = imag_x[7] - imag_x[9];
        /* RED */       real_v[3] += real_x[9];
        /* RED */       imag_v[3] += imag_x[9];
        /* RED */       real_v[10] = real_x[8] - real_x[9];
        /* RED */       imag_v[10] = imag_x[8] - imag_x[9];
        /* RED */       real_v[4] = real_x[10];
        /* RED */       imag_v[4] = imag_x[10];
        /* RED */       real_v[4] += real_x[11];
        /* RED */       imag_v[4] += imag_x[11];
        /* RED */       real_v[11] = real_x[10] - real_x[12];
        /* RED */       imag_v[11] = imag_x[10] - imag_x[12];
        /* RED */       real_v[4] += real_x[12];
        /* RED */       imag_v[4] += imag_x[12];
        /* RED */       real_v[12] = real_x[11] - real_x[12];
        /* RED */       imag_v[12] = imag_x[11] - imag_x[12];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULre  */       real_y[2] = real_x[2] * u[1];
        /* MULre  */       imag_y[2] = imag_x[2] * u[1];
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[3];
        /* D2  */       imag_v[0] = imag_x[3];
        /* D2  */       real_v[1] = real_x[4];
        /* D2  */       imag_v[1] = imag_x[4];
        /* D2  */       real_v[2] = real_x[3] + real_x[4];
        /* D2  */       imag_v[2] = imag_x[3] + imag_x[4];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[2];
        /* IMAG  */       imag_t = real_v[0] * u[2];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[3];
        /* IMAG  */       imag_t = real_v[1] * u[3];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[4];
        /* IMAG  */       imag_t = real_v[2] * u[4];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[3] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[3] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[4] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[4] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[5];
        /* D2  */       imag_v[0] = imag_x[5];
        /* D2  */       real_v[1] = real_x[6];
        /* D2  */       imag_v[1] = imag_x[6];
        /* D2  */       real_v[2] = real_x[5] + real_x[6];
        /* D2  */       imag_v[2] = imag_x[5] + imag_x[6];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[5];
        /* REAL  */       imag_v[0] = imag_v[0] * u[5];
        /* REAL  */       real_v[1] = real_v[1] * u[6];
        /* REAL  */       imag_v[1] = imag_v[1] * u[6];
        /* REAL  */       real_v[2] = real_v[2] * u[7];
        /* REAL  */       imag_v[2] = imag_v[2] * u[7];
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[5] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[5] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[6] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[6] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[7];
        /* D2  */       imag_v[0] = imag_x[7];
        /* D2  */       real_v[1] = real_x[8];
        /* D2  */       imag_v[1] = imag_x[8];
        /* D2  */       real_v[2] = real_x[7] + real_x[8];
        /* D2  */       imag_v[2] = imag_x[7] + imag_x[8];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[8];
        /* REAL  */       imag_v[0] = imag_v[0] * u[8];
        /* REAL  */       real_v[1] = real_v[1] * u[9];
        /* REAL  */       imag_v[1] = imag_v[1] * u[9];
        /* REAL  */       real_v[2] = real_v[2] * u[10];
        /* REAL  */       imag_v[2] = imag_v[2] * u[10];
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[7] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[7] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[8] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[8] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 2 */
        /* D2  */       real_v1[0] = real_x[9];
        /* D2  */       imag_v1[0] = imag_x[9];
        /* D2  */       real_v1[2] = real_x[11];
        /* D2  */       imag_v1[2] = imag_x[11];
        /* D2  */       real_v1[4] = real_x[9] + real_x[11];
        /* D2  */       imag_v1[4] = imag_x[9] + imag_x[11];
        /* D2  */       real_v1[1] = real_x[10];
        /* D2  */       imag_v1[1] = imag_x[10];
        /* D2  */       real_v1[3] = real_x[12];
        /* D2  */       imag_v1[3] = imag_x[12];
        /* D2  */       real_v1[5] = real_x[10] + real_x[12];
        /* D2  */       imag_v1[5] = imag_x[10] + imag_x[12];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 1 */
        /* D2  */       real_v[0] = real_v1[0];
        /* D2  */       imag_v[0] = imag_v1[0];
        /* D2  */       real_v[1] = real_v1[1];
        /* D2  */       imag_v[1] = imag_v1[1];
        /* D2  */       real_v[2] = real_v1[0] + real_v1[1];
        /* D2  */       imag_v[2] = imag_v1[0] + imag_v1[1];
        /* D2  */       real_v[3] = real_v1[2];
        /* D2  */       imag_v[3] = imag_v1[2];
        /* D2  */       real_v[4] = real_v1[3];
        /* D2  */       imag_v[4] = imag_v1[3];
        /* D2  */       real_v[5] = real_v1[2] + real_v1[3];
        /* D2  */       imag_v[5] = imag_v1[2] + imag_v1[3];
        /* D2  */       real_v[6] = real_v1[4];
        /* D2  */       imag_v[6] = imag_v1[4];
        /* D2  */       real_v[7] = real_v1[5];
        /* D2  */       imag_v[7] = imag_v1[5];
        /* D2  */       real_v[8] = real_v1[4] + real_v1[5];
        /* D2  */       imag_v[8] = imag_v1[4] + imag_v1[5];
        /* ID2I */     /* Exit */
     /* MARKER */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[11];
        /* IMAG  */       imag_t = real_v[0] * u[11];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[12];
        /* IMAG  */       imag_t = real_v[1] * u[12];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[13];
        /* IMAG  */       imag_t = real_v[2] * u[13];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[3] * u[14];
        /* IMAG  */       imag_t = real_v[3] * u[14];
        /* IMAG  */       real_v[3] = real_t;
        /* IMAG  */       imag_v[3] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[4] * u[15];
        /* IMAG  */       imag_t = real_v[4] * u[15];
        /* IMAG  */       real_v[4] = real_t;
        /* IMAG  */       imag_v[4] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[5] * u[16];
        /* IMAG  */       imag_t = real_v[5] * u[16];
        /* IMAG  */       real_v[5] = real_t;
        /* IMAG  */       imag_v[5] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[6] * u[17];
        /* IMAG  */       imag_t = real_v[6] * u[17];
        /* IMAG  */       real_v[6] = real_t;
        /* IMAG  */       imag_v[6] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[7] * u[18];
        /* IMAG  */       imag_t = real_v[7] * u[18];
        /* IMAG  */       real_v[7] = real_t;
        /* IMAG  */       imag_v[7] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[8] * u[19];
        /* IMAG  */       imag_t = real_v[8] * u[19];
        /* IMAG  */       real_v[8] = real_t;
        /* IMAG  */       imag_v[8] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 3  */
        /* D2t  */       real_v1[0] = real_v[0] + real_v[6];
        /* D2t  */       imag_v1[0] = imag_v[0] + imag_v[6];
        /* D2t  */       real_v1[3] = real_v[3] + real_v[6];
        /* D2t  */       imag_v1[3] = imag_v[3] + imag_v[6];
        /* D2t  */       real_v1[1] = real_v[1] + real_v[7];
        /* D2t  */       imag_v1[1] = imag_v[1] + imag_v[7];
        /* D2t  */       real_v1[4] = real_v[4] + real_v[7];
        /* D2t  */       imag_v1[4] = imag_v[4] + imag_v[7];
        /* D2t  */       real_v1[2] = real_v[2] + real_v[8];
        /* D2t  */       imag_v1[2] = imag_v[2] + imag_v[8];
        /* D2t  */       real_v1[5] = real_v[5] + real_v[8];
        /* D2t  */       imag_v1[5] = imag_v[5] + imag_v[8];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 1  */
        /* D2t  */       real_y[9] = real_v1[0] + real_v1[2];
        /* D2t  */       imag_y[9] = imag_v1[0] + imag_v1[2];
        /* D2t  */       real_y[10] = real_v1[1] + real_v1[2];
        /* D2t  */       imag_y[10] = imag_v1[1] + imag_v1[2];
        /* D2t  */       real_y[11] = real_v1[3] + real_v1[5];
        /* D2t  */       imag_y[11] = imag_v1[3] + imag_v1[5];
        /* D2t  */       real_y[12] = real_v1[4] + real_v1[5];
        /* D2t  */       imag_y[12] = imag_v1[4] + imag_v1[5];
        /* ID2It */     /* Exit */
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* tKRED */       /* tKRED entry */
          /* tRED */       /*  tRED entry p 3 a 4 c 1  */
        /* tRED */       real_v[3] = real_y[1];
        /* tRED */       imag_v[3] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[5];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[5];
        /* tRED */       real_v[3] -= real_y[5];
        /* tRED */       imag_v[3] -= imag_y[5];
        /* tRED */       real_v[2] = real_y[1] + real_y[6];
        /* tRED */       imag_v[2] = imag_y[1] + imag_y[6];
        /* tRED */       real_v[3] -= real_y[6];
        /* tRED */       imag_v[3] -= imag_y[6];
        /* tRED */       real_v[6] = real_y[2];
        /* tRED */       imag_v[6] = imag_y[2];
        /* tRED */       real_v[4] = real_y[2] + real_y[7];
        /* tRED */       imag_v[4] = imag_y[2] + imag_y[7];
        /* tRED */       real_v[6] -= real_y[7];
        /* tRED */       imag_v[6] -= imag_y[7];
        /* tRED */       real_v[5] = real_y[2] + real_y[8];
        /* tRED */       imag_v[5] = imag_y[2] + imag_y[8];
        /* tRED */       real_v[6] -= real_y[8];
        /* tRED */       imag_v[6] -= imag_y[8];
        /* tRED */       real_v[9] = real_y[3];
        /* tRED */       imag_v[9] = imag_y[3];
        /* tRED */       real_v[7] = real_y[3] + real_y[9];
        /* tRED */       imag_v[7] = imag_y[3] + imag_y[9];
        /* tRED */       real_v[9] -= real_y[9];
        /* tRED */       imag_v[9] -= imag_y[9];
        /* tRED */       real_v[8] = real_y[3] + real_y[10];
        /* tRED */       imag_v[8] = imag_y[3] + imag_y[10];
        /* tRED */       real_v[9] -= real_y[10];
        /* tRED */       imag_v[9] -= imag_y[10];
        /* tRED */       real_v[12] = real_y[4];
        /* tRED */       imag_v[12] = imag_y[4];
        /* tRED */       real_v[10] = real_y[4] + real_y[11];
        /* tRED */       imag_v[10] = imag_y[4] + imag_y[11];
        /* tRED */       real_v[12] -= real_y[11];
        /* tRED */       imag_v[12] -= imag_y[11];
        /* tRED */       real_v[11] = real_y[4] + real_y[12];
        /* tRED */       imag_v[11] = imag_y[4] + imag_y[12];
        /* tRED */       real_v[12] -= real_y[12];
        /* tRED */       imag_v[12] -= imag_y[12];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 3  */
        /* tRED */       real_v[4] = real_y[1];
        /* tRED */       imag_v[4] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[4];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[4];
        /* tRED */       real_v[4] -= real_y[4];
        /* tRED */       imag_v[4] -= imag_y[4];
        /* tRED */       real_v[5] = real_y[2];
        /* tRED */       imag_v[5] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[5];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[5];
        /* tRED */       real_v[5] -= real_y[5];
        /* tRED */       imag_v[5] -= imag_y[5];
        /* tRED */       real_v[6] = real_y[3];
        /* tRED */       imag_v[6] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[6];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[6];
        /* tRED */       real_v[6] -= real_y[6];
        /* tRED */       imag_v[6] -= imag_y[6];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 6  */
        /* tRED */       real_v[7] = real_y[1];
        /* tRED */       imag_v[7] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[7];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[7];
        /* tRED */       real_v[7] -= real_y[7];
        /* tRED */       imag_v[7] -= imag_y[7];
        /* tRED */       real_v[8] = real_y[2];
        /* tRED */       imag_v[8] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[8];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[8];
        /* tRED */       real_v[8] -= real_y[8];
        /* tRED */       imag_v[8] -= imag_y[8];
        /* tRED */       real_v[9] = real_y[3];
        /* tRED */       imag_v[9] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[9];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[9];
        /* tRED */       real_v[9] -= real_y[9];
        /* tRED */       imag_v[9] -= imag_y[9];
        /* tRED */       real_v[10] = real_y[4];
        /* tRED */       imag_v[10] = imag_y[4];
        /* tRED */       real_v[4] = real_y[4] + real_y[10];
        /* tRED */       imag_v[4] = imag_y[4] + imag_y[10];
        /* tRED */       real_v[10] -= real_y[10];
        /* tRED */       imag_v[10] -= imag_y[10];
        /* tRED */       real_v[11] = real_y[5];
        /* tRED */       imag_v[11] = imag_y[5];
        /* tRED */       real_v[5] = real_y[5] + real_y[11];
        /* tRED */       imag_v[5] = imag_y[5] + imag_y[11];
        /* tRED */       real_v[11] -= real_y[11];
        /* tRED */       imag_v[11] -= imag_y[11];
        /* tRED */       real_v[12] = real_y[6];
        /* tRED */       imag_v[12] = imag_y[6];
        /* tRED */       real_v[6] = real_y[6] + real_y[12];
        /* tRED */       imag_v[6] = imag_y[6] + imag_y[12];
        /* tRED */       real_v[12] -= real_y[12];
        /* tRED */       imag_v[12] -= imag_y[12];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       /*  tRED exit  */
      /* tKRED */       /*  tKRED exit */



        //data[op]		
        for (int px = 0; px < FFTLENGTH; px++)
        { 
            real_x[px] = real_y[active_op[px]]; 
            imag_x[px] = imag_y[active_op[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++)
        {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }
        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}


#undef FFTLENGTH
#define FFTLENGTH 17


void DFT17::Evaluate(Data *real, Data *imag)
{
    std::vector<s64> ind = indices;

    Data real_x[FFTLENGTH];
    Data imag_x[FFTLENGTH];
    Data real_v[27];
    Data imag_v[27];
    Data real_v1[18];
    Data imag_v1[18];
    Data real_y[FFTLENGTH];
    Data imag_y[FFTLENGTH];
    Data real_t, imag_t;

    for (int i = 0; i < count; i++)
    {

        for (int px = 0; px < FFTLENGTH; px++) {
            real_y[px] = real[ind[px]];
            imag_y[px] = imag[ind[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[ip[px]];
            imag_x[px] = imag_y[ip[px]];
        }


        //
    // DFT length 17
    //
  /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 8  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[9];
        /* RED */       imag_v[1] += imag_x[9];
        /* RED */       real_v[9] = real_x[1] - real_x[9];
        /* RED */       imag_v[9] = imag_x[1] - imag_x[9];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[10];
        /* RED */       imag_v[2] += imag_x[10];
        /* RED */       real_v[10] = real_x[2] - real_x[10];
        /* RED */       imag_v[10] = imag_x[2] - imag_x[10];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[11];
        /* RED */       imag_v[3] += imag_x[11];
        /* RED */       real_v[11] = real_x[3] - real_x[11];
        /* RED */       imag_v[11] = imag_x[3] - imag_x[11];
        /* RED */       real_v[4] = real_x[4];
        /* RED */       imag_v[4] = imag_x[4];
        /* RED */       real_v[4] += real_x[12];
        /* RED */       imag_v[4] += imag_x[12];
        /* RED */       real_v[12] = real_x[4] - real_x[12];
        /* RED */       imag_v[12] = imag_x[4] - imag_x[12];
        /* RED */       real_v[5] = real_x[5];
        /* RED */       imag_v[5] = imag_x[5];
        /* RED */       real_v[5] += real_x[13];
        /* RED */       imag_v[5] += imag_x[13];
        /* RED */       real_v[13] = real_x[5] - real_x[13];
        /* RED */       imag_v[13] = imag_x[5] - imag_x[13];
        /* RED */       real_v[6] = real_x[6];
        /* RED */       imag_v[6] = imag_x[6];
        /* RED */       real_v[6] += real_x[14];
        /* RED */       imag_v[6] += imag_x[14];
        /* RED */       real_v[14] = real_x[6] - real_x[14];
        /* RED */       imag_v[14] = imag_x[6] - imag_x[14];
        /* RED */       real_v[7] = real_x[7];
        /* RED */       imag_v[7] = imag_x[7];
        /* RED */       real_v[7] += real_x[15];
        /* RED */       imag_v[7] += imag_x[15];
        /* RED */       real_v[15] = real_x[7] - real_x[15];
        /* RED */       imag_v[15] = imag_x[7] - imag_x[15];
        /* RED */       real_v[8] = real_x[8];
        /* RED */       imag_v[8] = imag_x[8];
        /* RED */       real_v[8] += real_x[16];
        /* RED */       imag_v[8] += imag_x[16];
        /* RED */       real_v[16] = real_x[8] - real_x[16];
        /* RED */       imag_v[16] = imag_x[8] - imag_x[16];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       real_x[13] = real_v[13];
        /* RED */       imag_x[13] = imag_v[13];
        /* RED */       real_x[14] = real_v[14];
        /* RED */       imag_x[14] = imag_v[14];
        /* RED */       real_x[15] = real_v[15];
        /* RED */       imag_x[15] = imag_v[15];
        /* RED */       real_x[16] = real_v[16];
        /* RED */       imag_x[16] = imag_v[16];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 2  a 1 c 4  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[5];
        /* RED */       imag_v[1] += imag_x[5];
        /* RED */       real_v[5] = real_x[1] - real_x[5];
        /* RED */       imag_v[5] = imag_x[1] - imag_x[5];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[6];
        /* RED */       imag_v[2] += imag_x[6];
        /* RED */       real_v[6] = real_x[2] - real_x[6];
        /* RED */       imag_v[6] = imag_x[2] - imag_x[6];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[7];
        /* RED */       imag_v[3] += imag_x[7];
        /* RED */       real_v[7] = real_x[3] - real_x[7];
        /* RED */       imag_v[7] = imag_x[3] - imag_x[7];
        /* RED */       real_v[4] = real_x[4];
        /* RED */       imag_v[4] = imag_x[4];
        /* RED */       real_v[4] += real_x[8];
        /* RED */       imag_v[4] += imag_x[8];
        /* RED */       real_v[8] = real_x[4] - real_x[8];
        /* RED */       imag_v[8] = imag_x[4] - imag_x[8];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 2  a 1 c 2  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[3];
        /* RED */       imag_v[1] += imag_x[3];
        /* RED */       real_v[3] = real_x[1] - real_x[3];
        /* RED */       imag_v[3] = imag_x[1] - imag_x[3];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[4];
        /* RED */       imag_v[2] += imag_x[4];
        /* RED */       real_v[4] = real_x[2] - real_x[4];
        /* RED */       imag_v[4] = imag_x[2] - imag_x[4];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 2  a 1 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[2] = real_x[1] - real_x[2];
        /* RED */       imag_v[2] = imag_x[1] - imag_x[2];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULre  */       real_y[2] = real_x[2] * u[1];
        /* MULre  */       imag_y[2] = imag_x[2] * u[1];
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[3];
        /* D2  */       imag_v[0] = imag_x[3];
        /* D2  */       real_v[1] = real_x[4];
        /* D2  */       imag_v[1] = imag_x[4];
        /* D2  */       real_v[2] = real_x[3] + real_x[4];
        /* D2  */       imag_v[2] = imag_x[3] + imag_x[4];
        /* ID2I */     /* Exit */
     /* MARKER */
        /* REAL  */       real_v[0] = real_v[0] * u[2];
        /* REAL  */       imag_v[0] = imag_v[0] * u[2];
        /* REAL  */       real_v[1] = real_v[1] * u[3];
        /* REAL  */       imag_v[1] = imag_v[1] * u[3];
        /* REAL  */       real_v[2] = real_v[2] * u[4];
        /* REAL  */       imag_v[2] = imag_v[2] * u[4];
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[3] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[3] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[4] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[4] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 2 */
        /* D2  */       real_v[0] = real_x[5];
        /* D2  */       imag_v[0] = imag_x[5];
        /* D2  */       real_v[2] = real_x[7];
        /* D2  */       imag_v[2] = imag_x[7];
        /* D2  */       real_v[4] = real_x[5] + real_x[7];
        /* D2  */       imag_v[4] = imag_x[5] + imag_x[7];
        /* D2  */       real_v[1] = real_x[6];
        /* D2  */       imag_v[1] = imag_x[6];
        /* D2  */       real_v[3] = real_x[8];
        /* D2  */       imag_v[3] = imag_x[8];
        /* D2  */       real_v[5] = real_x[6] + real_x[8];
        /* D2  */       imag_v[5] = imag_x[6] + imag_x[8];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 1 */
        /* D2  */       real_v1[0] = real_v[0];
        /* D2  */       imag_v1[0] = imag_v[0];
        /* D2  */       real_v1[1] = real_v[1];
        /* D2  */       imag_v1[1] = imag_v[1];
        /* D2  */       real_v1[2] = real_v[0] + real_v[1];
        /* D2  */       imag_v1[2] = imag_v[0] + imag_v[1];
        /* D2  */       real_v1[3] = real_v[2];
        /* D2  */       imag_v1[3] = imag_v[2];
        /* D2  */       real_v1[4] = real_v[3];
        /* D2  */       imag_v1[4] = imag_v[3];
        /* D2  */       real_v1[5] = real_v[2] + real_v[3];
        /* D2  */       imag_v1[5] = imag_v[2] + imag_v[3];
        /* D2  */       real_v1[6] = real_v[4];
        /* D2  */       imag_v1[6] = imag_v[4];
        /* D2  */       real_v1[7] = real_v[5];
        /* D2  */       imag_v1[7] = imag_v[5];
        /* D2  */       real_v1[8] = real_v[4] + real_v[5];
        /* D2  */       imag_v1[8] = imag_v[4] + imag_v[5];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v1[0] = real_v1[0] * u[5];
        /* REAL  */       imag_v1[0] = imag_v1[0] * u[5];
        /* REAL  */       real_v1[1] = real_v1[1] * u[6];
        /* REAL  */       imag_v1[1] = imag_v1[1] * u[6];
        /* REAL  */       real_v1[2] = real_v1[2] * u[7];
        /* REAL  */       imag_v1[2] = imag_v1[2] * u[7];
        /* REAL  */       real_v1[3] = real_v1[3] * u[8];
        /* REAL  */       imag_v1[3] = imag_v1[3] * u[8];
        /* REAL  */       real_v1[4] = real_v1[4] * u[9];
        /* REAL  */       imag_v1[4] = imag_v1[4] * u[9];
        /* REAL  */       real_v1[5] = real_v1[5] * u[10];
        /* REAL  */       imag_v1[5] = imag_v1[5] * u[10];
        /* REAL  */       real_v1[6] = real_v1[6] * u[11];
        /* REAL  */       imag_v1[6] = imag_v1[6] * u[11];
        /* REAL  */       real_v1[7] = real_v1[7] * u[12];
        /* REAL  */       imag_v1[7] = imag_v1[7] * u[12];
        /* REAL  */       real_v1[8] = real_v1[8] * u[13];
        /* REAL  */       imag_v1[8] = imag_v1[8] * u[13];
        /* ID2It */     /* Entry m = 1  n = 3  */
        /* D2t  */       real_v[0] = real_v1[0] + real_v1[6];
        /* D2t  */       imag_v[0] = imag_v1[0] + imag_v1[6];
        /* D2t  */       real_v[3] = real_v1[3] + real_v1[6];
        /* D2t  */       imag_v[3] = imag_v1[3] + imag_v1[6];
        /* D2t  */       real_v[1] = real_v1[1] + real_v1[7];
        /* D2t  */       imag_v[1] = imag_v1[1] + imag_v1[7];
        /* D2t  */       real_v[4] = real_v1[4] + real_v1[7];
        /* D2t  */       imag_v[4] = imag_v1[4] + imag_v1[7];
        /* D2t  */       real_v[2] = real_v1[2] + real_v1[8];
        /* D2t  */       imag_v[2] = imag_v1[2] + imag_v1[8];
        /* D2t  */       real_v[5] = real_v1[5] + real_v1[8];
        /* D2t  */       imag_v[5] = imag_v1[5] + imag_v1[8];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 1  */
        /* D2t  */       real_y[5] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[5] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[6] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[6] = imag_v[1] + imag_v[2];
        /* D2t  */       real_y[7] = real_v[3] + real_v[5];
        /* D2t  */       imag_y[7] = imag_v[3] + imag_v[5];
        /* D2t  */       real_y[8] = real_v[4] + real_v[5];
        /* D2t  */       imag_y[8] = imag_v[4] + imag_v[5];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 4 */
        /* D2  */       real_v[0] = real_x[9];
        /* D2  */       imag_v[0] = imag_x[9];
        /* D2  */       real_v[4] = real_x[13];
        /* D2  */       imag_v[4] = imag_x[13];
        /* D2  */       real_v[8] = real_x[9] + real_x[13];
        /* D2  */       imag_v[8] = imag_x[9] + imag_x[13];
        /* D2  */       real_v[1] = real_x[10];
        /* D2  */       imag_v[1] = imag_x[10];
        /* D2  */       real_v[5] = real_x[14];
        /* D2  */       imag_v[5] = imag_x[14];
        /* D2  */       real_v[9] = real_x[10] + real_x[14];
        /* D2  */       imag_v[9] = imag_x[10] + imag_x[14];
        /* D2  */       real_v[2] = real_x[11];
        /* D2  */       imag_v[2] = imag_x[11];
        /* D2  */       real_v[6] = real_x[15];
        /* D2  */       imag_v[6] = imag_x[15];
        /* D2  */       real_v[10] = real_x[11] + real_x[15];
        /* D2  */       imag_v[10] = imag_x[11] + imag_x[15];
        /* D2  */       real_v[3] = real_x[12];
        /* D2  */       imag_v[3] = imag_x[12];
        /* D2  */       real_v[7] = real_x[16];
        /* D2  */       imag_v[7] = imag_x[16];
        /* D2  */       real_v[11] = real_x[12] + real_x[16];
        /* D2  */       imag_v[11] = imag_x[12] + imag_x[16];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 2 */
        /* D2  */       real_v1[0] = real_v[0];
        /* D2  */       imag_v1[0] = imag_v[0];
        /* D2  */       real_v1[2] = real_v[2];
        /* D2  */       imag_v1[2] = imag_v[2];
        /* D2  */       real_v1[4] = real_v[0] + real_v[2];
        /* D2  */       imag_v1[4] = imag_v[0] + imag_v[2];
        /* D2  */       real_v1[1] = real_v[1];
        /* D2  */       imag_v1[1] = imag_v[1];
        /* D2  */       real_v1[3] = real_v[3];
        /* D2  */       imag_v1[3] = imag_v[3];
        /* D2  */       real_v1[5] = real_v[1] + real_v[3];
        /* D2  */       imag_v1[5] = imag_v[1] + imag_v[3];
        /* D2  */       real_v1[6] = real_v[4];
        /* D2  */       imag_v1[6] = imag_v[4];
        /* D2  */       real_v1[8] = real_v[6];
        /* D2  */       imag_v1[8] = imag_v[6];
        /* D2  */       real_v1[10] = real_v[4] + real_v[6];
        /* D2  */       imag_v1[10] = imag_v[4] + imag_v[6];
        /* D2  */       real_v1[7] = real_v[5];
        /* D2  */       imag_v1[7] = imag_v[5];
        /* D2  */       real_v1[9] = real_v[7];
        /* D2  */       imag_v1[9] = imag_v[7];
        /* D2  */       real_v1[11] = real_v[5] + real_v[7];
        /* D2  */       imag_v1[11] = imag_v[5] + imag_v[7];
        /* D2  */       real_v1[12] = real_v[8];
        /* D2  */       imag_v1[12] = imag_v[8];
        /* D2  */       real_v1[14] = real_v[10];
        /* D2  */       imag_v1[14] = imag_v[10];
        /* D2  */       real_v1[16] = real_v[8] + real_v[10];
        /* D2  */       imag_v1[16] = imag_v[8] + imag_v[10];
        /* D2  */       real_v1[13] = real_v[9];
        /* D2  */       imag_v1[13] = imag_v[9];
        /* D2  */       real_v1[15] = real_v[11];
        /* D2  */       imag_v1[15] = imag_v[11];
        /* D2  */       real_v1[17] = real_v[9] + real_v[11];
        /* D2  */       imag_v1[17] = imag_v[9] + imag_v[11];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 9  n = 1 */
        /* D2  */       real_v[0] = real_v1[0];
        /* D2  */       imag_v[0] = imag_v1[0];
        /* D2  */       real_v[1] = real_v1[1];
        /* D2  */       imag_v[1] = imag_v1[1];
        /* D2  */       real_v[2] = real_v1[0] + real_v1[1];
        /* D2  */       imag_v[2] = imag_v1[0] + imag_v1[1];
        /* D2  */       real_v[3] = real_v1[2];
        /* D2  */       imag_v[3] = imag_v1[2];
        /* D2  */       real_v[4] = real_v1[3];
        /* D2  */       imag_v[4] = imag_v1[3];
        /* D2  */       real_v[5] = real_v1[2] + real_v1[3];
        /* D2  */       imag_v[5] = imag_v1[2] + imag_v1[3];
        /* D2  */       real_v[6] = real_v1[4];
        /* D2  */       imag_v[6] = imag_v1[4];
        /* D2  */       real_v[7] = real_v1[5];
        /* D2  */       imag_v[7] = imag_v1[5];
        /* D2  */       real_v[8] = real_v1[4] + real_v1[5];
        /* D2  */       imag_v[8] = imag_v1[4] + imag_v1[5];
        /* D2  */       real_v[9] = real_v1[6];
        /* D2  */       imag_v[9] = imag_v1[6];
        /* D2  */       real_v[10] = real_v1[7];
        /* D2  */       imag_v[10] = imag_v1[7];
        /* D2  */       real_v[11] = real_v1[6] + real_v1[7];
        /* D2  */       imag_v[11] = imag_v1[6] + imag_v1[7];
        /* D2  */       real_v[12] = real_v1[8];
        /* D2  */       imag_v[12] = imag_v1[8];
        /* D2  */       real_v[13] = real_v1[9];
        /* D2  */       imag_v[13] = imag_v1[9];
        /* D2  */       real_v[14] = real_v1[8] + real_v1[9];
        /* D2  */       imag_v[14] = imag_v1[8] + imag_v1[9];
        /* D2  */       real_v[15] = real_v1[10];
        /* D2  */       imag_v[15] = imag_v1[10];
        /* D2  */       real_v[16] = real_v1[11];
        /* D2  */       imag_v[16] = imag_v1[11];
        /* D2  */       real_v[17] = real_v1[10] + real_v1[11];
        /* D2  */       imag_v[17] = imag_v1[10] + imag_v1[11];
        /* D2  */       real_v[18] = real_v1[12];
        /* D2  */       imag_v[18] = imag_v1[12];
        /* D2  */       real_v[19] = real_v1[13];
        /* D2  */       imag_v[19] = imag_v1[13];
        /* D2  */       real_v[20] = real_v1[12] + real_v1[13];
        /* D2  */       imag_v[20] = imag_v1[12] + imag_v1[13];
        /* D2  */       real_v[21] = real_v1[14];
        /* D2  */       imag_v[21] = imag_v1[14];
        /* D2  */       real_v[22] = real_v1[15];
        /* D2  */       imag_v[22] = imag_v1[15];
        /* D2  */       real_v[23] = real_v1[14] + real_v1[15];
        /* D2  */       imag_v[23] = imag_v1[14] + imag_v1[15];
        /* D2  */       real_v[24] = real_v1[16];
        /* D2  */       imag_v[24] = imag_v1[16];
        /* D2  */       real_v[25] = real_v1[17];
        /* D2  */       imag_v[25] = imag_v1[17];
        /* D2  */       real_v[26] = real_v1[16] + real_v1[17];
        /* D2  */       imag_v[26] = imag_v1[16] + imag_v1[17];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[14];
        /* IMAG  */       imag_t = real_v[0] * u[14];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[15];
        /* IMAG  */       imag_t = real_v[1] * u[15];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[16];
        /* IMAG  */       imag_t = real_v[2] * u[16];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[3] * u[17];
        /* IMAG  */       imag_t = real_v[3] * u[17];
        /* IMAG  */       real_v[3] = real_t;
        /* IMAG  */       imag_v[3] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[4] * u[18];
        /* IMAG  */       imag_t = real_v[4] * u[18];
        /* IMAG  */       real_v[4] = real_t;
        /* IMAG  */       imag_v[4] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[5] * u[19];
        /* IMAG  */       imag_t = real_v[5] * u[19];
        /* IMAG  */       real_v[5] = real_t;
        /* IMAG  */       imag_v[5] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[6] * u[20];
        /* IMAG  */       imag_t = real_v[6] * u[20];
        /* IMAG  */       real_v[6] = real_t;
        /* IMAG  */       imag_v[6] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[7] * u[21];
        /* IMAG  */       imag_t = real_v[7] * u[21];
        /* IMAG  */       real_v[7] = real_t;
        /* IMAG  */       imag_v[7] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[8] * u[22];
        /* IMAG  */       imag_t = real_v[8] * u[22];
        /* IMAG  */       real_v[8] = real_t;
        /* IMAG  */       imag_v[8] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[9] * u[23];
        /* IMAG  */       imag_t = real_v[9] * u[23];
        /* IMAG  */       real_v[9] = real_t;
        /* IMAG  */       imag_v[9] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[10] * u[24];
        /* IMAG  */       imag_t = real_v[10] * u[24];
        /* IMAG  */       real_v[10] = real_t;
        /* IMAG  */       imag_v[10] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[11] * u[25];
        /* IMAG  */       imag_t = real_v[11] * u[25];
        /* IMAG  */       real_v[11] = real_t;
        /* IMAG  */       imag_v[11] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[12] * u[26];
        /* IMAG  */       imag_t = real_v[12] * u[26];
        /* IMAG  */       real_v[12] = real_t;
        /* IMAG  */       imag_v[12] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[13] * u[27];
        /* IMAG  */       imag_t = real_v[13] * u[27];
        /* IMAG  */       real_v[13] = real_t;
        /* IMAG  */       imag_v[13] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[14] * u[28];
        /* IMAG  */       imag_t = real_v[14] * u[28];
        /* IMAG  */       real_v[14] = real_t;
        /* IMAG  */       imag_v[14] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[15] * u[29];
        /* IMAG  */       imag_t = real_v[15] * u[29];
        /* IMAG  */       real_v[15] = real_t;
        /* IMAG  */       imag_v[15] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[16] * u[30];
        /* IMAG  */       imag_t = real_v[16] * u[30];
        /* IMAG  */       real_v[16] = real_t;
        /* IMAG  */       imag_v[16] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[17] * u[31];
        /* IMAG  */       imag_t = real_v[17] * u[31];
        /* IMAG  */       real_v[17] = real_t;
        /* IMAG  */       imag_v[17] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[18] * u[32];
        /* IMAG  */       imag_t = real_v[18] * u[32];
        /* IMAG  */       real_v[18] = real_t;
        /* IMAG  */       imag_v[18] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[19] * u[33];
        /* IMAG  */       imag_t = real_v[19] * u[33];
        /* IMAG  */       real_v[19] = real_t;
        /* IMAG  */       imag_v[19] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[20] * u[34];
        /* IMAG  */       imag_t = real_v[20] * u[34];
        /* IMAG  */       real_v[20] = real_t;
        /* IMAG  */       imag_v[20] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[21] * u[35];
        /* IMAG  */       imag_t = real_v[21] * u[35];
        /* IMAG  */       real_v[21] = real_t;
        /* IMAG  */       imag_v[21] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[22] * u[36];
        /* IMAG  */       imag_t = real_v[22] * u[36];
        /* IMAG  */       real_v[22] = real_t;
        /* IMAG  */       imag_v[22] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[23] * u[37];
        /* IMAG  */       imag_t = real_v[23] * u[37];
        /* IMAG  */       real_v[23] = real_t;
        /* IMAG  */       imag_v[23] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[24] * u[38];
        /* IMAG  */       imag_t = real_v[24] * u[38];
        /* IMAG  */       real_v[24] = real_t;
        /* IMAG  */       imag_v[24] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[25] * u[39];
        /* IMAG  */       imag_t = real_v[25] * u[39];
        /* IMAG  */       real_v[25] = real_t;
        /* IMAG  */       imag_v[25] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[26] * u[40];
        /* IMAG  */       imag_t = real_v[26] * u[40];
        /* IMAG  */       real_v[26] = real_t;
        /* IMAG  */       imag_v[26] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 9  */
        /* D2t  */       real_v1[0] = real_v[0] + real_v[18];
        /* D2t  */       imag_v1[0] = imag_v[0] + imag_v[18];
        /* D2t  */       real_v1[9] = real_v[9] + real_v[18];
        /* D2t  */       imag_v1[9] = imag_v[9] + imag_v[18];
        /* D2t  */       real_v1[1] = real_v[1] + real_v[19];
        /* D2t  */       imag_v1[1] = imag_v[1] + imag_v[19];
        /* D2t  */       real_v1[10] = real_v[10] + real_v[19];
        /* D2t  */       imag_v1[10] = imag_v[10] + imag_v[19];
        /* D2t  */       real_v1[2] = real_v[2] + real_v[20];
        /* D2t  */       imag_v1[2] = imag_v[2] + imag_v[20];
        /* D2t  */       real_v1[11] = real_v[11] + real_v[20];
        /* D2t  */       imag_v1[11] = imag_v[11] + imag_v[20];
        /* D2t  */       real_v1[3] = real_v[3] + real_v[21];
        /* D2t  */       imag_v1[3] = imag_v[3] + imag_v[21];
        /* D2t  */       real_v1[12] = real_v[12] + real_v[21];
        /* D2t  */       imag_v1[12] = imag_v[12] + imag_v[21];
        /* D2t  */       real_v1[4] = real_v[4] + real_v[22];
        /* D2t  */       imag_v1[4] = imag_v[4] + imag_v[22];
        /* D2t  */       real_v1[13] = real_v[13] + real_v[22];
        /* D2t  */       imag_v1[13] = imag_v[13] + imag_v[22];
        /* D2t  */       real_v1[5] = real_v[5] + real_v[23];
        /* D2t  */       imag_v1[5] = imag_v[5] + imag_v[23];
        /* D2t  */       real_v1[14] = real_v[14] + real_v[23];
        /* D2t  */       imag_v1[14] = imag_v[14] + imag_v[23];
        /* D2t  */       real_v1[6] = real_v[6] + real_v[24];
        /* D2t  */       imag_v1[6] = imag_v[6] + imag_v[24];
        /* D2t  */       real_v1[15] = real_v[15] + real_v[24];
        /* D2t  */       imag_v1[15] = imag_v[15] + imag_v[24];
        /* D2t  */       real_v1[7] = real_v[7] + real_v[25];
        /* D2t  */       imag_v1[7] = imag_v[7] + imag_v[25];
        /* D2t  */       real_v1[16] = real_v[16] + real_v[25];
        /* D2t  */       imag_v1[16] = imag_v[16] + imag_v[25];
        /* D2t  */       real_v1[8] = real_v[8] + real_v[26];
        /* D2t  */       imag_v1[8] = imag_v[8] + imag_v[26];
        /* D2t  */       real_v1[17] = real_v[17] + real_v[26];
        /* D2t  */       imag_v1[17] = imag_v[17] + imag_v[26];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 3  */
        /* D2t  */       real_v[0] = real_v1[0] + real_v1[6];
        /* D2t  */       imag_v[0] = imag_v1[0] + imag_v1[6];
        /* D2t  */       real_v[3] = real_v1[3] + real_v1[6];
        /* D2t  */       imag_v[3] = imag_v1[3] + imag_v1[6];
        /* D2t  */       real_v[1] = real_v1[1] + real_v1[7];
        /* D2t  */       imag_v[1] = imag_v1[1] + imag_v1[7];
        /* D2t  */       real_v[4] = real_v1[4] + real_v1[7];
        /* D2t  */       imag_v[4] = imag_v1[4] + imag_v1[7];
        /* D2t  */       real_v[2] = real_v1[2] + real_v1[8];
        /* D2t  */       imag_v[2] = imag_v1[2] + imag_v1[8];
        /* D2t  */       real_v[5] = real_v1[5] + real_v1[8];
        /* D2t  */       imag_v[5] = imag_v1[5] + imag_v1[8];
        /* D2t  */       real_v[6] = real_v1[9] + real_v1[15];
        /* D2t  */       imag_v[6] = imag_v1[9] + imag_v1[15];
        /* D2t  */       real_v[9] = real_v1[12] + real_v1[15];
        /* D2t  */       imag_v[9] = imag_v1[12] + imag_v1[15];
        /* D2t  */       real_v[7] = real_v1[10] + real_v1[16];
        /* D2t  */       imag_v[7] = imag_v1[10] + imag_v1[16];
        /* D2t  */       real_v[10] = real_v1[13] + real_v1[16];
        /* D2t  */       imag_v[10] = imag_v1[13] + imag_v1[16];
        /* D2t  */       real_v[8] = real_v1[11] + real_v1[17];
        /* D2t  */       imag_v[8] = imag_v1[11] + imag_v1[17];
        /* D2t  */       real_v[11] = real_v1[14] + real_v1[17];
        /* D2t  */       imag_v[11] = imag_v1[14] + imag_v1[17];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 4  n = 1  */
        /* D2t  */       real_y[9] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[9] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[10] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[10] = imag_v[1] + imag_v[2];
        /* D2t  */       real_y[11] = real_v[3] + real_v[5];
        /* D2t  */       imag_y[11] = imag_v[3] + imag_v[5];
        /* D2t  */       real_y[12] = real_v[4] + real_v[5];
        /* D2t  */       imag_y[12] = imag_v[4] + imag_v[5];
        /* D2t  */       real_y[13] = real_v[6] + real_v[8];
        /* D2t  */       imag_y[13] = imag_v[6] + imag_v[8];
        /* D2t  */       real_y[14] = real_v[7] + real_v[8];
        /* D2t  */       imag_y[14] = imag_v[7] + imag_v[8];
        /* D2t  */       real_y[15] = real_v[9] + real_v[11];
        /* D2t  */       imag_y[15] = imag_v[9] + imag_v[11];
        /* D2t  */       real_y[16] = real_v[10] + real_v[11];
        /* D2t  */       imag_y[16] = imag_v[10] + imag_v[11];
        /* ID2It */     /* Exit */
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* tKRED */       /* tKRED entry */
          /* tRED */       /*  tRED entry p 2 a 1 c 1  */
        /* tRED */       real_v[2] = real_y[1];
        /* tRED */       imag_v[2] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[2];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[2];
        /* tRED */       real_v[2] -= real_y[2];
        /* tRED */       imag_v[2] -= imag_y[2];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 2  */
        /* tRED */       real_v[3] = real_y[1];
        /* tRED */       imag_v[3] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[3];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[3];
        /* tRED */       real_v[3] -= real_y[3];
        /* tRED */       imag_v[3] -= imag_y[3];
        /* tRED */       real_v[4] = real_y[2];
        /* tRED */       imag_v[4] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[4];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[4];
        /* tRED */       real_v[4] -= real_y[4];
        /* tRED */       imag_v[4] -= imag_y[4];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 4  */
        /* tRED */       real_v[5] = real_y[1];
        /* tRED */       imag_v[5] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[5];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[5];
        /* tRED */       real_v[5] -= real_y[5];
        /* tRED */       imag_v[5] -= imag_y[5];
        /* tRED */       real_v[6] = real_y[2];
        /* tRED */       imag_v[6] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[6];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[6];
        /* tRED */       real_v[6] -= real_y[6];
        /* tRED */       imag_v[6] -= imag_y[6];
        /* tRED */       real_v[7] = real_y[3];
        /* tRED */       imag_v[7] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[7];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[7];
        /* tRED */       real_v[7] -= real_y[7];
        /* tRED */       imag_v[7] -= imag_y[7];
        /* tRED */       real_v[8] = real_y[4];
        /* tRED */       imag_v[8] = imag_y[4];
        /* tRED */       real_v[4] = real_y[4] + real_y[8];
        /* tRED */       imag_v[4] = imag_y[4] + imag_y[8];
        /* tRED */       real_v[8] -= real_y[8];
        /* tRED */       imag_v[8] -= imag_y[8];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 8  */
        /* tRED */       real_v[9] = real_y[1];
        /* tRED */       imag_v[9] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[9];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[9];
        /* tRED */       real_v[9] -= real_y[9];
        /* tRED */       imag_v[9] -= imag_y[9];
        /* tRED */       real_v[10] = real_y[2];
        /* tRED */       imag_v[10] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[10];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[10];
        /* tRED */       real_v[10] -= real_y[10];
        /* tRED */       imag_v[10] -= imag_y[10];
        /* tRED */       real_v[11] = real_y[3];
        /* tRED */       imag_v[11] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[11];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[11];
        /* tRED */       real_v[11] -= real_y[11];
        /* tRED */       imag_v[11] -= imag_y[11];
        /* tRED */       real_v[12] = real_y[4];
        /* tRED */       imag_v[12] = imag_y[4];
        /* tRED */       real_v[4] = real_y[4] + real_y[12];
        /* tRED */       imag_v[4] = imag_y[4] + imag_y[12];
        /* tRED */       real_v[12] -= real_y[12];
        /* tRED */       imag_v[12] -= imag_y[12];
        /* tRED */       real_v[13] = real_y[5];
        /* tRED */       imag_v[13] = imag_y[5];
        /* tRED */       real_v[5] = real_y[5] + real_y[13];
        /* tRED */       imag_v[5] = imag_y[5] + imag_y[13];
        /* tRED */       real_v[13] -= real_y[13];
        /* tRED */       imag_v[13] -= imag_y[13];
        /* tRED */       real_v[14] = real_y[6];
        /* tRED */       imag_v[14] = imag_y[6];
        /* tRED */       real_v[6] = real_y[6] + real_y[14];
        /* tRED */       imag_v[6] = imag_y[6] + imag_y[14];
        /* tRED */       real_v[14] -= real_y[14];
        /* tRED */       imag_v[14] -= imag_y[14];
        /* tRED */       real_v[15] = real_y[7];
        /* tRED */       imag_v[15] = imag_y[7];
        /* tRED */       real_v[7] = real_y[7] + real_y[15];
        /* tRED */       imag_v[7] = imag_y[7] + imag_y[15];
        /* tRED */       real_v[15] -= real_y[15];
        /* tRED */       imag_v[15] -= imag_y[15];
        /* tRED */       real_v[16] = real_y[8];
        /* tRED */       imag_v[16] = imag_y[8];
        /* tRED */       real_v[8] = real_y[8] + real_y[16];
        /* tRED */       imag_v[8] = imag_y[8] + imag_y[16];
        /* tRED */       real_v[16] -= real_y[16];
        /* tRED */       imag_v[16] -= imag_y[16];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       real_y[13] = real_v[13];
        /* tRED */       imag_y[13] = imag_v[13];
        /* tRED */       real_y[14] = real_v[14];
        /* tRED */       imag_y[14] = imag_v[14];
        /* tRED */       real_y[15] = real_v[15];
        /* tRED */       imag_y[15] = imag_v[15];
        /* tRED */       real_y[16] = real_v[16];
        /* tRED */       imag_y[16] = imag_v[16];
        /* tRED */       /*  tRED exit  */
      /* tKRED */       /*  tKRED exit */



    //data[op]		
        for (int px = 0; px < FFTLENGTH; px++)
        {
            real_x[px] = real_y[active_op[px]];
            imag_x[px] = imag_y[active_op[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++)
        {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }


        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}



#undef FFTLENGTH
#define FFTLENGTH 19

void DFT19::Evaluate(Data* real, Data* imag)
{
    std::vector<s64> ind = indices;


    Data real_x[FFTLENGTH];
    Data imag_x[FFTLENGTH];
    Data real_v[FFTLENGTH];
    Data imag_v[FFTLENGTH];
    Data real_v1[18];
    Data imag_v1[18];
    Data real_y[FFTLENGTH];
    Data imag_y[FFTLENGTH];
    Data real_t, imag_t;

    for (int i = 0; i < count; i++)
    {

        for (int px = 0; px < FFTLENGTH; px++) {
            real_y[px] = real[ind[px]];
            imag_y[px] = imag[ind[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[ip[px]];
            imag_x[px] = imag_y[ip[px]];
        }


        //
        // DFT length 19
        //
  /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 9  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[10];
        /* RED */       imag_v[1] += imag_x[10];
        /* RED */       real_v[10] = real_x[1] - real_x[10];
        /* RED */       imag_v[10] = imag_x[1] - imag_x[10];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[11];
        /* RED */       imag_v[2] += imag_x[11];
        /* RED */       real_v[11] = real_x[2] - real_x[11];
        /* RED */       imag_v[11] = imag_x[2] - imag_x[11];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[12];
        /* RED */       imag_v[3] += imag_x[12];
        /* RED */       real_v[12] = real_x[3] - real_x[12];
        /* RED */       imag_v[12] = imag_x[3] - imag_x[12];
        /* RED */       real_v[4] = real_x[4];
        /* RED */       imag_v[4] = imag_x[4];
        /* RED */       real_v[4] += real_x[13];
        /* RED */       imag_v[4] += imag_x[13];
        /* RED */       real_v[13] = real_x[4] - real_x[13];
        /* RED */       imag_v[13] = imag_x[4] - imag_x[13];
        /* RED */       real_v[5] = real_x[5];
        /* RED */       imag_v[5] = imag_x[5];
        /* RED */       real_v[5] += real_x[14];
        /* RED */       imag_v[5] += imag_x[14];
        /* RED */       real_v[14] = real_x[5] - real_x[14];
        /* RED */       imag_v[14] = imag_x[5] - imag_x[14];
        /* RED */       real_v[6] = real_x[6];
        /* RED */       imag_v[6] = imag_x[6];
        /* RED */       real_v[6] += real_x[15];
        /* RED */       imag_v[6] += imag_x[15];
        /* RED */       real_v[15] = real_x[6] - real_x[15];
        /* RED */       imag_v[15] = imag_x[6] - imag_x[15];
        /* RED */       real_v[7] = real_x[7];
        /* RED */       imag_v[7] = imag_x[7];
        /* RED */       real_v[7] += real_x[16];
        /* RED */       imag_v[7] += imag_x[16];
        /* RED */       real_v[16] = real_x[7] - real_x[16];
        /* RED */       imag_v[16] = imag_x[7] - imag_x[16];
        /* RED */       real_v[8] = real_x[8];
        /* RED */       imag_v[8] = imag_x[8];
        /* RED */       real_v[8] += real_x[17];
        /* RED */       imag_v[8] += imag_x[17];
        /* RED */       real_v[17] = real_x[8] - real_x[17];
        /* RED */       imag_v[17] = imag_x[8] - imag_x[17];
        /* RED */       real_v[9] = real_x[9];
        /* RED */       imag_v[9] = imag_x[9];
        /* RED */       real_v[9] += real_x[18];
        /* RED */       imag_v[9] += imag_x[18];
        /* RED */       real_v[18] = real_x[9] - real_x[18];
        /* RED */       imag_v[18] = imag_x[9] - imag_x[18];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       real_x[13] = real_v[13];
        /* RED */       imag_x[13] = imag_v[13];
        /* RED */       real_x[14] = real_v[14];
        /* RED */       imag_x[14] = imag_v[14];
        /* RED */       real_x[15] = real_v[15];
        /* RED */       imag_x[15] = imag_v[15];
        /* RED */       real_x[16] = real_v[16];
        /* RED */       imag_x[16] = imag_v[16];
        /* RED */       real_x[17] = real_v[17];
        /* RED */       imag_x[17] = imag_v[17];
        /* RED */       real_x[18] = real_v[18];
        /* RED */       imag_x[18] = imag_v[18];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 3  a 2 c 3  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[4];
        /* RED */       imag_v[1] += imag_x[4];
        /* RED */       real_v[7] = real_x[1] - real_x[7];
        /* RED */       imag_v[7] = imag_x[1] - imag_x[7];
        /* RED */       real_v[1] += real_x[7];
        /* RED */       imag_v[1] += imag_x[7];
        /* RED */       real_v[10] = real_x[4] - real_x[7];
        /* RED */       imag_v[10] = imag_x[4] - imag_x[7];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[5];
        /* RED */       imag_v[2] += imag_x[5];
        /* RED */       real_v[8] = real_x[2] - real_x[8];
        /* RED */       imag_v[8] = imag_x[2] - imag_x[8];
        /* RED */       real_v[2] += real_x[8];
        /* RED */       imag_v[2] += imag_x[8];
        /* RED */       real_v[11] = real_x[5] - real_x[8];
        /* RED */       imag_v[11] = imag_x[5] - imag_x[8];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[6];
        /* RED */       imag_v[3] += imag_x[6];
        /* RED */       real_v[9] = real_x[3] - real_x[9];
        /* RED */       imag_v[9] = imag_x[3] - imag_x[9];
        /* RED */       real_v[3] += real_x[9];
        /* RED */       imag_v[3] += imag_x[9];
        /* RED */       real_v[12] = real_x[6] - real_x[9];
        /* RED */       imag_v[12] = imag_x[6] - imag_x[9];
        /* RED */       real_v[4] = real_x[10];
        /* RED */       imag_v[4] = imag_x[10];
        /* RED */       real_v[4] += real_x[13];
        /* RED */       imag_v[4] += imag_x[13];
        /* RED */       real_v[13] = real_x[10] - real_x[16];
        /* RED */       imag_v[13] = imag_x[10] - imag_x[16];
        /* RED */       real_v[4] += real_x[16];
        /* RED */       imag_v[4] += imag_x[16];
        /* RED */       real_v[16] = real_x[13] - real_x[16];
        /* RED */       imag_v[16] = imag_x[13] - imag_x[16];
        /* RED */       real_v[5] = real_x[11];
        /* RED */       imag_v[5] = imag_x[11];
        /* RED */       real_v[5] += real_x[14];
        /* RED */       imag_v[5] += imag_x[14];
        /* RED */       real_v[14] = real_x[11] - real_x[17];
        /* RED */       imag_v[14] = imag_x[11] - imag_x[17];
        /* RED */       real_v[5] += real_x[17];
        /* RED */       imag_v[5] += imag_x[17];
        /* RED */       real_v[17] = real_x[14] - real_x[17];
        /* RED */       imag_v[17] = imag_x[14] - imag_x[17];
        /* RED */       real_v[6] = real_x[12];
        /* RED */       imag_v[6] = imag_x[12];
        /* RED */       real_v[6] += real_x[15];
        /* RED */       imag_v[6] += imag_x[15];
        /* RED */       real_v[15] = real_x[12] - real_x[18];
        /* RED */       imag_v[15] = imag_x[12] - imag_x[18];
        /* RED */       real_v[6] += real_x[18];
        /* RED */       imag_v[6] += imag_x[18];
        /* RED */       real_v[18] = real_x[15] - real_x[18];
        /* RED */       imag_v[18] = imag_x[15] - imag_x[18];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       real_x[13] = real_v[13];
        /* RED */       imag_x[13] = imag_v[13];
        /* RED */       real_x[14] = real_v[14];
        /* RED */       imag_x[14] = imag_v[14];
        /* RED */       real_x[15] = real_v[15];
        /* RED */       imag_x[15] = imag_v[15];
        /* RED */       real_x[16] = real_v[16];
        /* RED */       imag_x[16] = imag_v[16];
        /* RED */       real_x[17] = real_v[17];
        /* RED */       imag_x[17] = imag_v[17];
        /* RED */       real_x[18] = real_v[18];
        /* RED */       imag_x[18] = imag_v[18];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 3  a 2 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[3] = real_x[1] - real_x[3];
        /* RED */       imag_v[3] = imag_x[1] - imag_x[3];
        /* RED */       real_v[1] += real_x[3];
        /* RED */       imag_v[1] += imag_x[3];
        /* RED */       real_v[4] = real_x[2] - real_x[3];
        /* RED */       imag_v[4] = imag_x[2] - imag_x[3];
        /* RED */       real_v[2] = real_x[4];
        /* RED */       imag_v[2] = imag_x[4];
        /* RED */       real_v[2] += real_x[5];
        /* RED */       imag_v[2] += imag_x[5];
        /* RED */       real_v[5] = real_x[4] - real_x[6];
        /* RED */       imag_v[5] = imag_x[4] - imag_x[6];
        /* RED */       real_v[2] += real_x[6];
        /* RED */       imag_v[2] += imag_x[6];
        /* RED */       real_v[6] = real_x[5] - real_x[6];
        /* RED */       imag_v[6] = imag_x[5] - imag_x[6];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULim  */       real_y[2] = -1 * imag_x[2] * u[1];
        /* MULim  */       imag_y[2] = real_x[2] * u[1];
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[3];
        /* D2  */       imag_v[0] = imag_x[3];
        /* D2  */       real_v[1] = real_x[4];
        /* D2  */       imag_v[1] = imag_x[4];
        /* D2  */       real_v[2] = real_x[3] + real_x[4];
        /* D2  */       imag_v[2] = imag_x[3] + imag_x[4];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[2];
        /* REAL  */       imag_v[0] = imag_v[0] * u[2];
        /* REAL  */       real_v[1] = real_v[1] * u[3];
        /* REAL  */       imag_v[1] = imag_v[1] * u[3];
        /* REAL  */       real_v[2] = real_v[2] * u[4];
        /* REAL  */       imag_v[2] = imag_v[2] * u[4];
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[3] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[3] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[4] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[4] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[5];
        /* D2  */       imag_v[0] = imag_x[5];
        /* D2  */       real_v[1] = real_x[6];
        /* D2  */       imag_v[1] = imag_x[6];
        /* D2  */       real_v[2] = real_x[5] + real_x[6];
        /* D2  */       imag_v[2] = imag_x[5] + imag_x[6];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[5];
        /* IMAG  */       imag_t = real_v[0] * u[5];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[6];
        /* IMAG  */       imag_t = real_v[1] * u[6];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[7];
        /* IMAG  */       imag_t = real_v[2] * u[7];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[5] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[5] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[6] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[6] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID3I */     /* Entry m = 2  n = 1 */
        /* D3  */ { Data real_a, real_b, imag_a, imag_b;
        /* D3  */       real_a = real_x[8] + real_x[9];
        /* D3  */       imag_a = imag_x[8] + imag_x[9];
        /* D3  */       real_b = real_x[9] - real_x[8];
        /* D3  */       imag_b = imag_x[9] - imag_x[8];
        /* D3  */       real_v1[0] = real_x[7];
        /* D3  */       imag_v1[0] = imag_x[7];
        /* D3  */       real_v1[1] = real_x[7] + real_a;
        /* D3  */       imag_v1[1] = imag_x[7] + imag_a;
        /* D3  */       real_v1[2] = real_x[7] + real_b;
        /* D3  */       imag_v1[2] = imag_x[7] + imag_b;
        /* D3  */       real_v1[3] = real_a + real_a + real_b + real_v1[1];
        /* D3  */       imag_v1[3] = imag_a + imag_a + imag_b + imag_v1[1];
        /* D2  */       real_v1[4] = real_x[9];
        /* D2  */       imag_v1[4] = imag_x[9];
        /* D3  */       }
        /* D3  */ { Data real_a, real_b,imag_a, imag_b;
        /* D3  */       real_a = real_x[11] + real_x[12];
        /* D3  */       imag_a = imag_x[11] + imag_x[12];
        /* D3  */       real_b = real_x[12] - real_x[11];
        /* D3  */       imag_b = imag_x[12] - imag_x[11];
        /* D3  */       real_v1[5] = real_x[10];
        /* D3  */       imag_v1[5] = imag_x[10];
        /* D3  */       real_v1[6] = real_x[10] + real_a;
        /* D3  */       imag_v1[6] = imag_x[10] + imag_a;
        /* D3  */       real_v1[7] = real_x[10] + real_b;
        /* D3  */       imag_v1[7] = imag_x[10] + imag_b;
        /* D3  */       real_v1[8] = real_a + real_a + real_b + real_v1[6];
        /* D3  */       imag_v1[8] = imag_a + imag_a + imag_b + imag_v1[6];
        /* D2  */       real_v1[9] = real_x[12];
        /* D2  */       imag_v1[9] = imag_x[12];
        /* D3  */       }
        /* ID3I */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 5 */
        /* D2  */       real_v[0] = real_v1[0];
        /* D2  */       imag_v[0] = imag_v1[0];
        /* D2  */       real_v[5] = real_v1[5];
        /* D2  */       imag_v[5] = imag_v1[5];
        /* D2  */       real_v[10] = real_v1[0] + real_v1[5];
        /* D2  */       imag_v[10] = imag_v1[0] + imag_v1[5];
        /* D2  */       real_v[1] = real_v1[1];
        /* D2  */       imag_v[1] = imag_v1[1];
        /* D2  */       real_v[6] = real_v1[6];
        /* D2  */       imag_v[6] = imag_v1[6];
        /* D2  */       real_v[11] = real_v1[1] + real_v1[6];
        /* D2  */       imag_v[11] = imag_v1[1] + imag_v1[6];
        /* D2  */       real_v[2] = real_v1[2];
        /* D2  */       imag_v[2] = imag_v1[2];
        /* D2  */       real_v[7] = real_v1[7];
        /* D2  */       imag_v[7] = imag_v1[7];
        /* D2  */       real_v[12] = real_v1[2] + real_v1[7];
        /* D2  */       imag_v[12] = imag_v1[2] + imag_v1[7];
        /* D2  */       real_v[3] = real_v1[3];
        /* D2  */       imag_v[3] = imag_v1[3];
        /* D2  */       real_v[8] = real_v1[8];
        /* D2  */       imag_v[8] = imag_v1[8];
        /* D2  */       real_v[13] = real_v1[3] + real_v1[8];
        /* D2  */       imag_v[13] = imag_v1[3] + imag_v1[8];
        /* D2  */       real_v[4] = real_v1[4];
        /* D2  */       imag_v[4] = imag_v1[4];
        /* D2  */       real_v[9] = real_v1[9];
        /* D2  */       imag_v[9] = imag_v1[9];
        /* D2  */       real_v[14] = real_v1[4] + real_v1[9];
        /* D2  */       imag_v[14] = imag_v1[4] + imag_v1[9];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[8];
        /* REAL  */       imag_v[0] = imag_v[0] * u[8];
        /* REAL  */       real_v[1] = real_v[1] * u[9];
        /* REAL  */       imag_v[1] = imag_v[1] * u[9];
        /* REAL  */       real_v[2] = real_v[2] * u[10];
        /* REAL  */       imag_v[2] = imag_v[2] * u[10];
        /* REAL  */       real_v[3] = real_v[3] * u[11];
        /* REAL  */       imag_v[3] = imag_v[3] * u[11];
        /* REAL  */       real_v[4] = real_v[4] * u[12];
        /* REAL  */       imag_v[4] = imag_v[4] * u[12];
        /* REAL  */       real_v[5] = real_v[5] * u[13];
        /* REAL  */       imag_v[5] = imag_v[5] * u[13];
        /* REAL  */       real_v[6] = real_v[6] * u[14];
        /* REAL  */       imag_v[6] = imag_v[6] * u[14];
        /* REAL  */       real_v[7] = real_v[7] * u[15];
        /* REAL  */       imag_v[7] = imag_v[7] * u[15];
        /* REAL  */       real_v[8] = real_v[8] * u[16];
        /* REAL  */       imag_v[8] = imag_v[8] * u[16];
        /* REAL  */       real_v[9] = real_v[9] * u[17];
        /* REAL  */       imag_v[9] = imag_v[9] * u[17];
        /* REAL  */       real_v[10] = real_v[10] * u[18];
        /* REAL  */       imag_v[10] = imag_v[10] * u[18];
        /* REAL  */       real_v[11] = real_v[11] * u[19];
        /* REAL  */       imag_v[11] = imag_v[11] * u[19];
        /* REAL  */       real_v[12] = real_v[12] * u[20];
        /* REAL  */       imag_v[12] = imag_v[12] * u[20];
        /* REAL  */       real_v[13] = real_v[13] * u[21];
        /* REAL  */       imag_v[13] = imag_v[13] * u[21];
        /* REAL  */       real_v[14] = real_v[14] * u[22];
        /* REAL  */       imag_v[14] = imag_v[14] * u[22];
        /* ID2It */     /* Entry m = 1  n = 5  */
        /* D2t  */       real_v1[0] = real_v[0] + real_v[10];
        /* D2t  */       imag_v1[0] = imag_v[0] + imag_v[10];
        /* D2t  */       real_v1[5] = real_v[5] + real_v[10];
        /* D2t  */       imag_v1[5] = imag_v[5] + imag_v[10];
        /* D2t  */       real_v1[1] = real_v[1] + real_v[11];
        /* D2t  */       imag_v1[1] = imag_v[1] + imag_v[11];
        /* D2t  */       real_v1[6] = real_v[6] + real_v[11];
        /* D2t  */       imag_v1[6] = imag_v[6] + imag_v[11];
        /* D2t  */       real_v1[2] = real_v[2] + real_v[12];
        /* D2t  */       imag_v1[2] = imag_v[2] + imag_v[12];
        /* D2t  */       real_v1[7] = real_v[7] + real_v[12];
        /* D2t  */       imag_v1[7] = imag_v[7] + imag_v[12];
        /* D2t  */       real_v1[3] = real_v[3] + real_v[13];
        /* D2t  */       imag_v1[3] = imag_v[3] + imag_v[13];
        /* D2t  */       real_v1[8] = real_v[8] + real_v[13];
        /* D2t  */       imag_v1[8] = imag_v[8] + imag_v[13];
        /* D2t  */       real_v1[4] = real_v[4] + real_v[14];
        /* D2t  */       imag_v1[4] = imag_v[4] + imag_v[14];
        /* D2t  */       real_v1[9] = real_v[9] + real_v[14];
        /* D2t  */       imag_v1[9] = imag_v[9] + imag_v[14];
        /* ID2It */     /* Exit */
        /* ID3It */     /* Entry m = 2  n = 1  */
        /* D3t  */ { Data real_a, imag_a ;
        /* D3t  */       real_y[7] = real_v1[1] + real_v1[2] + real_v1[3];
        /* D3t  */       imag_y[7] = imag_v1[1] + imag_v1[2] + imag_v1[3];
        /* D3t  */       real_a = real_v1[3] + real_v1[3];
        /* D3t  */       imag_a = imag_v1[3] + imag_v1[3];
        /* D3t  */       real_y[8] = real_v1[1] - real_v1[2] + real_a;
        /* D3t  */       imag_y[8] = imag_v1[1] - imag_v1[2] + imag_a;
        /* D3t  */       real_y[9] = real_y[7] + real_v1[3] + real_a;
        /* D3t  */       imag_y[9] = imag_y[7] + imag_v1[3] + imag_a;
        /* D3t  */       real_y[7] = real_y[7] + real_v1[0];
        /* D3t  */       imag_y[7] = imag_y[7] + imag_v1[0];
        /* D3t  */       real_y[9] = real_y[9] + real_v1[4];
        /* D3t  */       imag_y[9] = imag_y[9] + imag_v1[4];
        /* D3t  */       }
        /* D3t  */ { Data real_a, imag_a ;
        /* D3t  */       real_y[10] = real_v1[6] + real_v1[7] + real_v1[8];
        /* D3t  */       imag_y[10] = imag_v1[6] + imag_v1[7] + imag_v1[8];
        /* D3t  */       real_a = real_v1[8] + real_v1[8];
        /* D3t  */       imag_a = imag_v1[8] + imag_v1[8];
        /* D3t  */       real_y[11] = real_v1[6] - real_v1[7] + real_a;
        /* D3t  */       imag_y[11] = imag_v1[6] - imag_v1[7] + imag_a;
        /* D3t  */       real_y[12] = real_y[10] + real_v1[8] + real_a;
        /* D3t  */       imag_y[12] = imag_y[10] + imag_v1[8] + imag_a;
        /* D3t  */       real_y[10] = real_y[10] + real_v1[5];
        /* D3t  */       imag_y[10] = imag_y[10] + imag_v1[5];
        /* D3t  */       real_y[12] = real_y[12] + real_v1[9];
        /* D3t  */       imag_y[12] = imag_y[12] + imag_v1[9];
        /* D3t  */       }
        /* ID3It */     /* Exit */
        /* ID3I */     /* Entry m = 2  n = 1 */
        /* D3  */ { Data real_a, real_b, imag_a, imag_b;
        /* D3  */       real_a = real_x[14] + real_x[15];
        /* D3  */       imag_a = imag_x[14] + imag_x[15];
        /* D3  */       real_b = real_x[15] - real_x[14];
        /* D3  */       imag_b = imag_x[15] - imag_x[14];
        /* D3  */       real_v[0] = real_x[13];
        /* D3  */       imag_v[0] = imag_x[13];
        /* D3  */       real_v[1] = real_x[13] + real_a;
        /* D3  */       imag_v[1] = imag_x[13] + imag_a;
        /* D3  */       real_v[2] = real_x[13] + real_b;
        /* D3  */       imag_v[2] = imag_x[13] + imag_b;
        /* D3  */       real_v[3] = real_a + real_a + real_b + real_v[1];
        /* D3  */       imag_v[3] = imag_a + imag_a + imag_b + imag_v[1];
        /* D2  */       real_v[4] = real_x[15];
        /* D2  */       imag_v[4] = imag_x[15];
        /* D3  */       }
        /* D3  */ { Data real_a, real_b,  imag_a, imag_b;
        /* D3  */       real_a = real_x[17] + real_x[18];
        /* D3  */       imag_a = imag_x[17] + imag_x[18];
        /* D3  */       real_b = real_x[18] - real_x[17];
        /* D3  */       imag_b = imag_x[18] - imag_x[17];
        /* D3  */       real_v[5] = real_x[16];
        /* D3  */       imag_v[5] = imag_x[16];
        /* D3  */       real_v[6] = real_x[16] + real_a;
        /* D3  */       imag_v[6] = imag_x[16] + imag_a;
        /* D3  */       real_v[7] = real_x[16] + real_b;
        /* D3  */       imag_v[7] = imag_x[16] + imag_b;
        /* D3  */       real_v[8] = real_a + real_a + real_b + real_v[6];
        /* D3  */       imag_v[8] = imag_a + imag_a + imag_b + imag_v[6];
        /* D2  */       real_v[9] = real_x[18];
        /* D2  */       imag_v[9] = imag_x[18];
        /* D3  */       }
        /* ID3I */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 5 */
        /* D2  */       real_v1[0] = real_v[0];
        /* D2  */       imag_v1[0] = imag_v[0];
        /* D2  */       real_v1[5] = real_v[5];
        /* D2  */       imag_v1[5] = imag_v[5];
        /* D2  */       real_v1[10] = real_v[0] + real_v[5];
        /* D2  */       imag_v1[10] = imag_v[0] + imag_v[5];
        /* D2  */       real_v1[1] = real_v[1];
        /* D2  */       imag_v1[1] = imag_v[1];
        /* D2  */       real_v1[6] = real_v[6];
        /* D2  */       imag_v1[6] = imag_v[6];
        /* D2  */       real_v1[11] = real_v[1] + real_v[6];
        /* D2  */       imag_v1[11] = imag_v[1] + imag_v[6];
        /* D2  */       real_v1[2] = real_v[2];
        /* D2  */       imag_v1[2] = imag_v[2];
        /* D2  */       real_v1[7] = real_v[7];
        /* D2  */       imag_v1[7] = imag_v[7];
        /* D2  */       real_v1[12] = real_v[2] + real_v[7];
        /* D2  */       imag_v1[12] = imag_v[2] + imag_v[7];
        /* D2  */       real_v1[3] = real_v[3];
        /* D2  */       imag_v1[3] = imag_v[3];
        /* D2  */       real_v1[8] = real_v[8];
        /* D2  */       imag_v1[8] = imag_v[8];
        /* D2  */       real_v1[13] = real_v[3] + real_v[8];
        /* D2  */       imag_v1[13] = imag_v[3] + imag_v[8];
        /* D2  */       real_v1[4] = real_v[4];
        /* D2  */       imag_v1[4] = imag_v[4];
        /* D2  */       real_v1[9] = real_v[9];
        /* D2  */       imag_v1[9] = imag_v[9];
        /* D2  */       real_v1[14] = real_v[4] + real_v[9];
        /* D2  */       imag_v1[14] = imag_v[4] + imag_v[9];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v1[0] * u[23];
        /* IMAG  */       imag_t = real_v1[0] * u[23];
        /* IMAG  */       real_v1[0] = real_t;
        /* IMAG  */       imag_v1[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[1] * u[24];
        /* IMAG  */       imag_t = real_v1[1] * u[24];
        /* IMAG  */       real_v1[1] = real_t;
        /* IMAG  */       imag_v1[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[2] * u[25];
        /* IMAG  */       imag_t = real_v1[2] * u[25];
        /* IMAG  */       real_v1[2] = real_t;
        /* IMAG  */       imag_v1[2] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[3] * u[26];
        /* IMAG  */       imag_t = real_v1[3] * u[26];
        /* IMAG  */       real_v1[3] = real_t;
        /* IMAG  */       imag_v1[3] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[4] * u[27];
        /* IMAG  */       imag_t = real_v1[4] * u[27];
        /* IMAG  */       real_v1[4] = real_t;
        /* IMAG  */       imag_v1[4] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[5] * u[28];
        /* IMAG  */       imag_t = real_v1[5] * u[28];
        /* IMAG  */       real_v1[5] = real_t;
        /* IMAG  */       imag_v1[5] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[6] * u[29];
        /* IMAG  */       imag_t = real_v1[6] * u[29];
        /* IMAG  */       real_v1[6] = real_t;
        /* IMAG  */       imag_v1[6] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[7] * u[30];
        /* IMAG  */       imag_t = real_v1[7] * u[30];
        /* IMAG  */       real_v1[7] = real_t;
        /* IMAG  */       imag_v1[7] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[8] * u[31];
        /* IMAG  */       imag_t = real_v1[8] * u[31];
        /* IMAG  */       real_v1[8] = real_t;
        /* IMAG  */       imag_v1[8] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[9] * u[32];
        /* IMAG  */       imag_t = real_v1[9] * u[32];
        /* IMAG  */       real_v1[9] = real_t;
        /* IMAG  */       imag_v1[9] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[10] * u[33];
        /* IMAG  */       imag_t = real_v1[10] * u[33];
        /* IMAG  */       real_v1[10] = real_t;
        /* IMAG  */       imag_v1[10] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[11] * u[34];
        /* IMAG  */       imag_t = real_v1[11] * u[34];
        /* IMAG  */       real_v1[11] = real_t;
        /* IMAG  */       imag_v1[11] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[12] * u[35];
        /* IMAG  */       imag_t = real_v1[12] * u[35];
        /* IMAG  */       real_v1[12] = real_t;
        /* IMAG  */       imag_v1[12] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[13] * u[36];
        /* IMAG  */       imag_t = real_v1[13] * u[36];
        /* IMAG  */       real_v1[13] = real_t;
        /* IMAG  */       imag_v1[13] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v1[14] * u[37];
        /* IMAG  */       imag_t = real_v1[14] * u[37];
        /* IMAG  */       real_v1[14] = real_t;
        /* IMAG  */       imag_v1[14] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 5  */
        /* D2t  */       real_v[0] = real_v1[0] + real_v1[10];
        /* D2t  */       imag_v[0] = imag_v1[0] + imag_v1[10];
        /* D2t  */       real_v[5] = real_v1[5] + real_v1[10];
        /* D2t  */       imag_v[5] = imag_v1[5] + imag_v1[10];
        /* D2t  */       real_v[1] = real_v1[1] + real_v1[11];
        /* D2t  */       imag_v[1] = imag_v1[1] + imag_v1[11];
        /* D2t  */       real_v[6] = real_v1[6] + real_v1[11];
        /* D2t  */       imag_v[6] = imag_v1[6] + imag_v1[11];
        /* D2t  */       real_v[2] = real_v1[2] + real_v1[12];
        /* D2t  */       imag_v[2] = imag_v1[2] + imag_v1[12];
        /* D2t  */       real_v[7] = real_v1[7] + real_v1[12];
        /* D2t  */       imag_v[7] = imag_v1[7] + imag_v1[12];
        /* D2t  */       real_v[3] = real_v1[3] + real_v1[13];
        /* D2t  */       imag_v[3] = imag_v1[3] + imag_v1[13];
        /* D2t  */       real_v[8] = real_v1[8] + real_v1[13];
        /* D2t  */       imag_v[8] = imag_v1[8] + imag_v1[13];
        /* D2t  */       real_v[4] = real_v1[4] + real_v1[14];
        /* D2t  */       imag_v[4] = imag_v1[4] + imag_v1[14];
        /* D2t  */       real_v[9] = real_v1[9] + real_v1[14];
        /* D2t  */       imag_v[9] = imag_v1[9] + imag_v1[14];
        /* ID2It */     /* Exit */
        /* ID3It */     /* Entry m = 2  n = 1  */
        /* D3t  */ { Data real_a, imag_a ;
        /* D3t  */       real_y[13] = real_v[1] + real_v[2] + real_v[3];
        /* D3t  */       imag_y[13] = imag_v[1] + imag_v[2] + imag_v[3];
        /* D3t  */       real_a = real_v[3] + real_v[3];
        /* D3t  */       imag_a = imag_v[3] + imag_v[3];
        /* D3t  */       real_y[14] = real_v[1] - real_v[2] + real_a;
        /* D3t  */       imag_y[14] = imag_v[1] - imag_v[2] + imag_a;
        /* D3t  */       real_y[15] = real_y[13] + real_v[3] + real_a;
        /* D3t  */       imag_y[15] = imag_y[13] + imag_v[3] + imag_a;
        /* D3t  */       real_y[13] = real_y[13] + real_v[0];
        /* D3t  */       imag_y[13] = imag_y[13] + imag_v[0];
        /* D3t  */       real_y[15] = real_y[15] + real_v[4];
        /* D3t  */       imag_y[15] = imag_y[15] + imag_v[4];
        /* D3t  */       }
        /* D3t  */ { Data real_a, imag_a ;
        /* D3t  */       real_y[16] = real_v[6] + real_v[7] + real_v[8];
        /* D3t  */       imag_y[16] = imag_v[6] + imag_v[7] + imag_v[8];
        /* D3t  */       real_a = real_v[8] + real_v[8];
        /* D3t  */       imag_a = imag_v[8] + imag_v[8];
        /* D3t  */       real_y[17] = real_v[6] - real_v[7] + real_a;
        /* D3t  */       imag_y[17] = imag_v[6] - imag_v[7] + imag_a;
        /* D3t  */       real_y[18] = real_y[16] + real_v[8] + real_a;
        /* D3t  */       imag_y[18] = imag_y[16] + imag_v[8] + imag_a;
        /* D3t  */       real_y[16] = real_y[16] + real_v[5];
        /* D3t  */       imag_y[16] = imag_y[16] + imag_v[5];
        /* D3t  */       real_y[18] = real_y[18] + real_v[9];
        /* D3t  */       imag_y[18] = imag_y[18] + imag_v[9];
        /* D3t  */       }
        /* ID3It */     /* Exit */
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* MARKER */
         /* tKRED */       /* tKRED entry */
           /* tRED */       /*  tRED entry p 3 a 2 c 1  */
        /* tRED */       real_v[3] = real_y[1];
        /* tRED */       imag_v[3] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[3];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[3];
        /* tRED */       real_v[3] -= real_y[3];
        /* tRED */       imag_v[3] -= imag_y[3];
        /* tRED */       real_v[2] = real_y[1] + real_y[4];
        /* tRED */       imag_v[2] = imag_y[1] + imag_y[4];
        /* tRED */       real_v[3] -= real_y[4];
        /* tRED */       imag_v[3] -= imag_y[4];
        /* tRED */       real_v[6] = real_y[2];
        /* tRED */       imag_v[6] = imag_y[2];
        /* tRED */       real_v[4] = real_y[2] + real_y[5];
        /* tRED */       imag_v[4] = imag_y[2] + imag_y[5];
        /* tRED */       real_v[6] -= real_y[5];
        /* tRED */       imag_v[6] -= imag_y[5];
        /* tRED */       real_v[5] = real_y[2] + real_y[6];
        /* tRED */       imag_v[5] = imag_y[2] + imag_y[6];
        /* tRED */       real_v[6] -= real_y[6];
        /* tRED */       imag_v[6] -= imag_y[6];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 3 a 2 c 3  */
        /* tRED */       real_v[7] = real_y[1];
        /* tRED */       imag_v[7] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[7];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[7];
        /* tRED */       real_v[7] -= real_y[7];
        /* tRED */       imag_v[7] -= imag_y[7];
        /* tRED */       real_v[4] = real_y[1] + real_y[10];
        /* tRED */       imag_v[4] = imag_y[1] + imag_y[10];
        /* tRED */       real_v[7] -= real_y[10];
        /* tRED */       imag_v[7] -= imag_y[10];
        /* tRED */       real_v[8] = real_y[2];
        /* tRED */       imag_v[8] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[8];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[8];
        /* tRED */       real_v[8] -= real_y[8];
        /* tRED */       imag_v[8] -= imag_y[8];
        /* tRED */       real_v[5] = real_y[2] + real_y[11];
        /* tRED */       imag_v[5] = imag_y[2] + imag_y[11];
        /* tRED */       real_v[8] -= real_y[11];
        /* tRED */       imag_v[8] -= imag_y[11];
        /* tRED */       real_v[9] = real_y[3];
        /* tRED */       imag_v[9] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[9];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[9];
        /* tRED */       real_v[9] -= real_y[9];
        /* tRED */       imag_v[9] -= imag_y[9];
        /* tRED */       real_v[6] = real_y[3] + real_y[12];
        /* tRED */       imag_v[6] = imag_y[3] + imag_y[12];
        /* tRED */       real_v[9] -= real_y[12];
        /* tRED */       imag_v[9] -= imag_y[12];
        /* tRED */       real_v[16] = real_y[4];
        /* tRED */       imag_v[16] = imag_y[4];
        /* tRED */       real_v[10] = real_y[4] + real_y[13];
        /* tRED */       imag_v[10] = imag_y[4] + imag_y[13];
        /* tRED */       real_v[16] -= real_y[13];
        /* tRED */       imag_v[16] -= imag_y[13];
        /* tRED */       real_v[13] = real_y[4] + real_y[16];
        /* tRED */       imag_v[13] = imag_y[4] + imag_y[16];
        /* tRED */       real_v[16] -= real_y[16];
        /* tRED */       imag_v[16] -= imag_y[16];
        /* tRED */       real_v[17] = real_y[5];
        /* tRED */       imag_v[17] = imag_y[5];
        /* tRED */       real_v[11] = real_y[5] + real_y[14];
        /* tRED */       imag_v[11] = imag_y[5] + imag_y[14];
        /* tRED */       real_v[17] -= real_y[14];
        /* tRED */       imag_v[17] -= imag_y[14];
        /* tRED */       real_v[14] = real_y[5] + real_y[17];
        /* tRED */       imag_v[14] = imag_y[5] + imag_y[17];
        /* tRED */       real_v[17] -= real_y[17];
        /* tRED */       imag_v[17] -= imag_y[17];
        /* tRED */       real_v[18] = real_y[6];
        /* tRED */       imag_v[18] = imag_y[6];
        /* tRED */       real_v[12] = real_y[6] + real_y[15];
        /* tRED */       imag_v[12] = imag_y[6] + imag_y[15];
        /* tRED */       real_v[18] -= real_y[15];
        /* tRED */       imag_v[18] -= imag_y[15];
        /* tRED */       real_v[15] = real_y[6] + real_y[18];
        /* tRED */       imag_v[15] = imag_y[6] + imag_y[18];
        /* tRED */       real_v[18] -= real_y[18];
        /* tRED */       imag_v[18] -= imag_y[18];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       real_y[13] = real_v[13];
        /* tRED */       imag_y[13] = imag_v[13];
        /* tRED */       real_y[14] = real_v[14];
        /* tRED */       imag_y[14] = imag_v[14];
        /* tRED */       real_y[15] = real_v[15];
        /* tRED */       imag_y[15] = imag_v[15];
        /* tRED */       real_y[16] = real_v[16];
        /* tRED */       imag_y[16] = imag_v[16];
        /* tRED */       real_y[17] = real_v[17];
        /* tRED */       imag_y[17] = imag_v[17];
        /* tRED */       real_y[18] = real_v[18];
        /* tRED */       imag_y[18] = imag_v[18];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 9  */
        /* tRED */       real_v[10] = real_y[1];
        /* tRED */       imag_v[10] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[10];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[10];
        /* tRED */       real_v[10] -= real_y[10];
        /* tRED */       imag_v[10] -= imag_y[10];
        /* tRED */       real_v[11] = real_y[2];
        /* tRED */       imag_v[11] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[11];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[11];
        /* tRED */       real_v[11] -= real_y[11];
        /* tRED */       imag_v[11] -= imag_y[11];
        /* tRED */       real_v[12] = real_y[3];
        /* tRED */       imag_v[12] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[12];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[12];
        /* tRED */       real_v[12] -= real_y[12];
        /* tRED */       imag_v[12] -= imag_y[12];
        /* tRED */       real_v[13] = real_y[4];
        /* tRED */       imag_v[13] = imag_y[4];
        /* tRED */       real_v[4] = real_y[4] + real_y[13];
        /* tRED */       imag_v[4] = imag_y[4] + imag_y[13];
        /* tRED */       real_v[13] -= real_y[13];
        /* tRED */       imag_v[13] -= imag_y[13];
        /* tRED */       real_v[14] = real_y[5];
        /* tRED */       imag_v[14] = imag_y[5];
        /* tRED */       real_v[5] = real_y[5] + real_y[14];
        /* tRED */       imag_v[5] = imag_y[5] + imag_y[14];
        /* tRED */       real_v[14] -= real_y[14];
        /* tRED */       imag_v[14] -= imag_y[14];
        /* tRED */       real_v[15] = real_y[6];
        /* tRED */       imag_v[15] = imag_y[6];
        /* tRED */       real_v[6] = real_y[6] + real_y[15];
        /* tRED */       imag_v[6] = imag_y[6] + imag_y[15];
        /* tRED */       real_v[15] -= real_y[15];
        /* tRED */       imag_v[15] -= imag_y[15];
        /* tRED */       real_v[16] = real_y[7];
        /* tRED */       imag_v[16] = imag_y[7];
        /* tRED */       real_v[7] = real_y[7] + real_y[16];
        /* tRED */       imag_v[7] = imag_y[7] + imag_y[16];
        /* tRED */       real_v[16] -= real_y[16];
        /* tRED */       imag_v[16] -= imag_y[16];
        /* tRED */       real_v[17] = real_y[8];
        /* tRED */       imag_v[17] = imag_y[8];
        /* tRED */       real_v[8] = real_y[8] + real_y[17];
        /* tRED */       imag_v[8] = imag_y[8] + imag_y[17];
        /* tRED */       real_v[17] -= real_y[17];
        /* tRED */       imag_v[17] -= imag_y[17];
        /* tRED */       real_v[18] = real_y[9];
        /* tRED */       imag_v[18] = imag_y[9];
        /* tRED */       real_v[9] = real_y[9] + real_y[18];
        /* tRED */       imag_v[9] = imag_y[9] + imag_y[18];
        /* tRED */       real_v[18] -= real_y[18];
        /* tRED */       imag_v[18] -= imag_y[18];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       real_y[13] = real_v[13];
        /* tRED */       imag_y[13] = imag_v[13];
        /* tRED */       real_y[14] = real_v[14];
        /* tRED */       imag_y[14] = imag_v[14];
        /* tRED */       real_y[15] = real_v[15];
        /* tRED */       imag_y[15] = imag_v[15];
        /* tRED */       real_y[16] = real_v[16];
        /* tRED */       imag_y[16] = imag_v[16];
        /* tRED */       real_y[17] = real_v[17];
        /* tRED */       imag_y[17] = imag_v[17];
        /* tRED */       real_y[18] = real_v[18];
        /* tRED */       imag_y[18] = imag_v[18];
        /* tRED */       /*  tRED exit  */
      /* tKRED */       /*  tKRED exit */

        //data[op]		
        for (int px = 0; px < FFTLENGTH; px++)
        {
            real_x[px] = real_y[active_op[px]];
            imag_x[px] = imag_y[active_op[px]];
        }


        for (int px = 0; px < FFTLENGTH; px++)
        {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }


        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}



#undef FFTLENGTH 
#define FFTLENGTH 31

void DFT31::Evaluate(Data  *real, Data *imag)
{
    std::vector<s64> ind = indices;


    Data real_x[FFTLENGTH];
    Data imag_x[FFTLENGTH];
    Data real_v[FFTLENGTH];
    Data imag_v[FFTLENGTH];
    Data real_v1[18];
    Data imag_v1[18];
    Data real_y[FFTLENGTH];
    Data imag_y[FFTLENGTH];
    Data real_t, imag_t;

    for (int i = 0; i < count; i++)
    {

        for (int px = 0; px < FFTLENGTH; px++) {
            real_y[px] = real[ind[px]];
            imag_y[px] = imag[ind[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++) {
            real_x[px] = real_y[ip[px]];
            imag_x[px] = imag_y[ip[px]];
        }


        //
        // DFT length 31
        //
  /* KRED */       /* KRED entry */
    /* RED */       /*  RED entry p 2  a 1 c 15  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[16];
        /* RED */       imag_v[1] += imag_x[16];
        /* RED */       real_v[16] = real_x[1] - real_x[16];
        /* RED */       imag_v[16] = imag_x[1] - imag_x[16];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[17];
        /* RED */       imag_v[2] += imag_x[17];
        /* RED */       real_v[17] = real_x[2] - real_x[17];
        /* RED */       imag_v[17] = imag_x[2] - imag_x[17];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[18];
        /* RED */       imag_v[3] += imag_x[18];
        /* RED */       real_v[18] = real_x[3] - real_x[18];
        /* RED */       imag_v[18] = imag_x[3] - imag_x[18];
        /* RED */       real_v[4] = real_x[4];
        /* RED */       imag_v[4] = imag_x[4];
        /* RED */       real_v[4] += real_x[19];
        /* RED */       imag_v[4] += imag_x[19];
        /* RED */       real_v[19] = real_x[4] - real_x[19];
        /* RED */       imag_v[19] = imag_x[4] - imag_x[19];
        /* RED */       real_v[5] = real_x[5];
        /* RED */       imag_v[5] = imag_x[5];
        /* RED */       real_v[5] += real_x[20];
        /* RED */       imag_v[5] += imag_x[20];
        /* RED */       real_v[20] = real_x[5] - real_x[20];
        /* RED */       imag_v[20] = imag_x[5] - imag_x[20];
        /* RED */       real_v[6] = real_x[6];
        /* RED */       imag_v[6] = imag_x[6];
        /* RED */       real_v[6] += real_x[21];
        /* RED */       imag_v[6] += imag_x[21];
        /* RED */       real_v[21] = real_x[6] - real_x[21];
        /* RED */       imag_v[21] = imag_x[6] - imag_x[21];
        /* RED */       real_v[7] = real_x[7];
        /* RED */       imag_v[7] = imag_x[7];
        /* RED */       real_v[7] += real_x[22];
        /* RED */       imag_v[7] += imag_x[22];
        /* RED */       real_v[22] = real_x[7] - real_x[22];
        /* RED */       imag_v[22] = imag_x[7] - imag_x[22];
        /* RED */       real_v[8] = real_x[8];
        /* RED */       imag_v[8] = imag_x[8];
        /* RED */       real_v[8] += real_x[23];
        /* RED */       imag_v[8] += imag_x[23];
        /* RED */       real_v[23] = real_x[8] - real_x[23];
        /* RED */       imag_v[23] = imag_x[8] - imag_x[23];
        /* RED */       real_v[9] = real_x[9];
        /* RED */       imag_v[9] = imag_x[9];
        /* RED */       real_v[9] += real_x[24];
        /* RED */       imag_v[9] += imag_x[24];
        /* RED */       real_v[24] = real_x[9] - real_x[24];
        /* RED */       imag_v[24] = imag_x[9] - imag_x[24];
        /* RED */       real_v[10] = real_x[10];
        /* RED */       imag_v[10] = imag_x[10];
        /* RED */       real_v[10] += real_x[25];
        /* RED */       imag_v[10] += imag_x[25];
        /* RED */       real_v[25] = real_x[10] - real_x[25];
        /* RED */       imag_v[25] = imag_x[10] - imag_x[25];
        /* RED */       real_v[11] = real_x[11];
        /* RED */       imag_v[11] = imag_x[11];
        /* RED */       real_v[11] += real_x[26];
        /* RED */       imag_v[11] += imag_x[26];
        /* RED */       real_v[26] = real_x[11] - real_x[26];
        /* RED */       imag_v[26] = imag_x[11] - imag_x[26];
        /* RED */       real_v[12] = real_x[12];
        /* RED */       imag_v[12] = imag_x[12];
        /* RED */       real_v[12] += real_x[27];
        /* RED */       imag_v[12] += imag_x[27];
        /* RED */       real_v[27] = real_x[12] - real_x[27];
        /* RED */       imag_v[27] = imag_x[12] - imag_x[27];
        /* RED */       real_v[13] = real_x[13];
        /* RED */       imag_v[13] = imag_x[13];
        /* RED */       real_v[13] += real_x[28];
        /* RED */       imag_v[13] += imag_x[28];
        /* RED */       real_v[28] = real_x[13] - real_x[28];
        /* RED */       imag_v[28] = imag_x[13] - imag_x[28];
        /* RED */       real_v[14] = real_x[14];
        /* RED */       imag_v[14] = imag_x[14];
        /* RED */       real_v[14] += real_x[29];
        /* RED */       imag_v[14] += imag_x[29];
        /* RED */       real_v[29] = real_x[14] - real_x[29];
        /* RED */       imag_v[29] = imag_x[14] - imag_x[29];
        /* RED */       real_v[15] = real_x[15];
        /* RED */       imag_v[15] = imag_x[15];
        /* RED */       real_v[15] += real_x[30];
        /* RED */       imag_v[15] += imag_x[30];
        /* RED */       real_v[30] = real_x[15] - real_x[30];
        /* RED */       imag_v[30] = imag_x[15] - imag_x[30];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       real_x[13] = real_v[13];
        /* RED */       imag_x[13] = imag_v[13];
        /* RED */       real_x[14] = real_v[14];
        /* RED */       imag_x[14] = imag_v[14];
        /* RED */       real_x[15] = real_v[15];
        /* RED */       imag_x[15] = imag_v[15];
        /* RED */       real_x[16] = real_v[16];
        /* RED */       imag_x[16] = imag_v[16];
        /* RED */       real_x[17] = real_v[17];
        /* RED */       imag_x[17] = imag_v[17];
        /* RED */       real_x[18] = real_v[18];
        /* RED */       imag_x[18] = imag_v[18];
        /* RED */       real_x[19] = real_v[19];
        /* RED */       imag_x[19] = imag_v[19];
        /* RED */       real_x[20] = real_v[20];
        /* RED */       imag_x[20] = imag_v[20];
        /* RED */       real_x[21] = real_v[21];
        /* RED */       imag_x[21] = imag_v[21];
        /* RED */       real_x[22] = real_v[22];
        /* RED */       imag_x[22] = imag_v[22];
        /* RED */       real_x[23] = real_v[23];
        /* RED */       imag_x[23] = imag_v[23];
        /* RED */       real_x[24] = real_v[24];
        /* RED */       imag_x[24] = imag_v[24];
        /* RED */       real_x[25] = real_v[25];
        /* RED */       imag_x[25] = imag_v[25];
        /* RED */       real_x[26] = real_v[26];
        /* RED */       imag_x[26] = imag_v[26];
        /* RED */       real_x[27] = real_v[27];
        /* RED */       imag_x[27] = imag_v[27];
        /* RED */       real_x[28] = real_v[28];
        /* RED */       imag_x[28] = imag_v[28];
        /* RED */       real_x[29] = real_v[29];
        /* RED */       imag_x[29] = imag_v[29];
        /* RED */       real_x[30] = real_v[30];
        /* RED */       imag_x[30] = imag_v[30];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 3  a 2 c 5  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[6];
        /* RED */       imag_v[1] += imag_x[6];
        /* RED */       real_v[11] = real_x[1] - real_x[11];
        /* RED */       imag_v[11] = imag_x[1] - imag_x[11];
        /* RED */       real_v[1] += real_x[11];
        /* RED */       imag_v[1] += imag_x[11];
        /* RED */       real_v[16] = real_x[6] - real_x[11];
        /* RED */       imag_v[16] = imag_x[6] - imag_x[11];
        /* RED */       real_v[2] = real_x[2];
        /* RED */       imag_v[2] = imag_x[2];
        /* RED */       real_v[2] += real_x[7];
        /* RED */       imag_v[2] += imag_x[7];
        /* RED */       real_v[12] = real_x[2] - real_x[12];
        /* RED */       imag_v[12] = imag_x[2] - imag_x[12];
        /* RED */       real_v[2] += real_x[12];
        /* RED */       imag_v[2] += imag_x[12];
        /* RED */       real_v[17] = real_x[7] - real_x[12];
        /* RED */       imag_v[17] = imag_x[7] - imag_x[12];
        /* RED */       real_v[3] = real_x[3];
        /* RED */       imag_v[3] = imag_x[3];
        /* RED */       real_v[3] += real_x[8];
        /* RED */       imag_v[3] += imag_x[8];
        /* RED */       real_v[13] = real_x[3] - real_x[13];
        /* RED */       imag_v[13] = imag_x[3] - imag_x[13];
        /* RED */       real_v[3] += real_x[13];
        /* RED */       imag_v[3] += imag_x[13];
        /* RED */       real_v[18] = real_x[8] - real_x[13];
        /* RED */       imag_v[18] = imag_x[8] - imag_x[13];
        /* RED */       real_v[4] = real_x[4];
        /* RED */       imag_v[4] = imag_x[4];
        /* RED */       real_v[4] += real_x[9];
        /* RED */       imag_v[4] += imag_x[9];
        /* RED */       real_v[14] = real_x[4] - real_x[14];
        /* RED */       imag_v[14] = imag_x[4] - imag_x[14];
        /* RED */       real_v[4] += real_x[14];
        /* RED */       imag_v[4] += imag_x[14];
        /* RED */       real_v[19] = real_x[9] - real_x[14];
        /* RED */       imag_v[19] = imag_x[9] - imag_x[14];
        /* RED */       real_v[5] = real_x[5];
        /* RED */       imag_v[5] = imag_x[5];
        /* RED */       real_v[5] += real_x[10];
        /* RED */       imag_v[5] += imag_x[10];
        /* RED */       real_v[15] = real_x[5] - real_x[15];
        /* RED */       imag_v[15] = imag_x[5] - imag_x[15];
        /* RED */       real_v[5] += real_x[15];
        /* RED */       imag_v[5] += imag_x[15];
        /* RED */       real_v[20] = real_x[10] - real_x[15];
        /* RED */       imag_v[20] = imag_x[10] - imag_x[15];
        /* RED */       real_v[6] = real_x[16];
        /* RED */       imag_v[6] = imag_x[16];
        /* RED */       real_v[6] += real_x[21];
        /* RED */       imag_v[6] += imag_x[21];
        /* RED */       real_v[21] = real_x[16] - real_x[26];
        /* RED */       imag_v[21] = imag_x[16] - imag_x[26];
        /* RED */       real_v[6] += real_x[26];
        /* RED */       imag_v[6] += imag_x[26];
        /* RED */       real_v[26] = real_x[21] - real_x[26];
        /* RED */       imag_v[26] = imag_x[21] - imag_x[26];
        /* RED */       real_v[7] = real_x[17];
        /* RED */       imag_v[7] = imag_x[17];
        /* RED */       real_v[7] += real_x[22];
        /* RED */       imag_v[7] += imag_x[22];
        /* RED */       real_v[22] = real_x[17] - real_x[27];
        /* RED */       imag_v[22] = imag_x[17] - imag_x[27];
        /* RED */       real_v[7] += real_x[27];
        /* RED */       imag_v[7] += imag_x[27];
        /* RED */       real_v[27] = real_x[22] - real_x[27];
        /* RED */       imag_v[27] = imag_x[22] - imag_x[27];
        /* RED */       real_v[8] = real_x[18];
        /* RED */       imag_v[8] = imag_x[18];
        /* RED */       real_v[8] += real_x[23];
        /* RED */       imag_v[8] += imag_x[23];
        /* RED */       real_v[23] = real_x[18] - real_x[28];
        /* RED */       imag_v[23] = imag_x[18] - imag_x[28];
        /* RED */       real_v[8] += real_x[28];
        /* RED */       imag_v[8] += imag_x[28];
        /* RED */       real_v[28] = real_x[23] - real_x[28];
        /* RED */       imag_v[28] = imag_x[23] - imag_x[28];
        /* RED */       real_v[9] = real_x[19];
        /* RED */       imag_v[9] = imag_x[19];
        /* RED */       real_v[9] += real_x[24];
        /* RED */       imag_v[9] += imag_x[24];
        /* RED */       real_v[24] = real_x[19] - real_x[29];
        /* RED */       imag_v[24] = imag_x[19] - imag_x[29];
        /* RED */       real_v[9] += real_x[29];
        /* RED */       imag_v[9] += imag_x[29];
        /* RED */       real_v[29] = real_x[24] - real_x[29];
        /* RED */       imag_v[29] = imag_x[24] - imag_x[29];
        /* RED */       real_v[10] = real_x[20];
        /* RED */       imag_v[10] = imag_x[20];
        /* RED */       real_v[10] += real_x[25];
        /* RED */       imag_v[10] += imag_x[25];
        /* RED */       real_v[25] = real_x[20] - real_x[30];
        /* RED */       imag_v[25] = imag_x[20] - imag_x[30];
        /* RED */       real_v[10] += real_x[30];
        /* RED */       imag_v[10] += imag_x[30];
        /* RED */       real_v[30] = real_x[25] - real_x[30];
        /* RED */       imag_v[30] = imag_x[25] - imag_x[30];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       real_x[13] = real_v[13];
        /* RED */       imag_x[13] = imag_v[13];
        /* RED */       real_x[14] = real_v[14];
        /* RED */       imag_x[14] = imag_v[14];
        /* RED */       real_x[15] = real_v[15];
        /* RED */       imag_x[15] = imag_v[15];
        /* RED */       real_x[16] = real_v[16];
        /* RED */       imag_x[16] = imag_v[16];
        /* RED */       real_x[17] = real_v[17];
        /* RED */       imag_x[17] = imag_v[17];
        /* RED */       real_x[18] = real_v[18];
        /* RED */       imag_x[18] = imag_v[18];
        /* RED */       real_x[19] = real_v[19];
        /* RED */       imag_x[19] = imag_v[19];
        /* RED */       real_x[20] = real_v[20];
        /* RED */       imag_x[20] = imag_v[20];
        /* RED */       real_x[21] = real_v[21];
        /* RED */       imag_x[21] = imag_v[21];
        /* RED */       real_x[22] = real_v[22];
        /* RED */       imag_x[22] = imag_v[22];
        /* RED */       real_x[23] = real_v[23];
        /* RED */       imag_x[23] = imag_v[23];
        /* RED */       real_x[24] = real_v[24];
        /* RED */       imag_x[24] = imag_v[24];
        /* RED */       real_x[25] = real_v[25];
        /* RED */       imag_x[25] = imag_v[25];
        /* RED */       real_x[26] = real_v[26];
        /* RED */       imag_x[26] = imag_v[26];
        /* RED */       real_x[27] = real_v[27];
        /* RED */       imag_x[27] = imag_v[27];
        /* RED */       real_x[28] = real_v[28];
        /* RED */       imag_x[28] = imag_v[28];
        /* RED */       real_x[29] = real_v[29];
        /* RED */       imag_x[29] = imag_v[29];
        /* RED */       real_x[30] = real_v[30];
        /* RED */       imag_x[30] = imag_v[30];
        /* RED */       /*  RED exit */
        /* RED */       /*  RED entry p 5  a 6 c 1  */
        /* RED */       real_v[1] = real_x[1];
        /* RED */       imag_v[1] = imag_x[1];
        /* RED */       real_v[1] += real_x[2];
        /* RED */       imag_v[1] += imag_x[2];
        /* RED */       real_v[7] = real_x[1] - real_x[5];
        /* RED */       imag_v[7] = imag_x[1] - imag_x[5];
        /* RED */       real_v[1] += real_x[3];
        /* RED */       imag_v[1] += imag_x[3];
        /* RED */       real_v[8] = real_x[2] - real_x[5];
        /* RED */       imag_v[8] = imag_x[2] - imag_x[5];
        /* RED */       real_v[1] += real_x[4];
        /* RED */       imag_v[1] += imag_x[4];
        /* RED */       real_v[9] = real_x[3] - real_x[5];
        /* RED */       imag_v[9] = imag_x[3] - imag_x[5];
        /* RED */       real_v[1] += real_x[5];
        /* RED */       imag_v[1] += imag_x[5];
        /* RED */       real_v[10] = real_x[4] - real_x[5];
        /* RED */       imag_v[10] = imag_x[4] - imag_x[5];
        /* RED */       real_v[2] = real_x[6];
        /* RED */       imag_v[2] = imag_x[6];
        /* RED */       real_v[2] += real_x[7];
        /* RED */       imag_v[2] += imag_x[7];
        /* RED */       real_v[11] = real_x[6] - real_x[10];
        /* RED */       imag_v[11] = imag_x[6] - imag_x[10];
        /* RED */       real_v[2] += real_x[8];
        /* RED */       imag_v[2] += imag_x[8];
        /* RED */       real_v[12] = real_x[7] - real_x[10];
        /* RED */       imag_v[12] = imag_x[7] - imag_x[10];
        /* RED */       real_v[2] += real_x[9];
        /* RED */       imag_v[2] += imag_x[9];
        /* RED */       real_v[13] = real_x[8] - real_x[10];
        /* RED */       imag_v[13] = imag_x[8] - imag_x[10];
        /* RED */       real_v[2] += real_x[10];
        /* RED */       imag_v[2] += imag_x[10];
        /* RED */       real_v[14] = real_x[9] - real_x[10];
        /* RED */       imag_v[14] = imag_x[9] - imag_x[10];
        /* RED */       real_v[3] = real_x[11];
        /* RED */       imag_v[3] = imag_x[11];
        /* RED */       real_v[3] += real_x[12];
        /* RED */       imag_v[3] += imag_x[12];
        /* RED */       real_v[15] = real_x[11] - real_x[15];
        /* RED */       imag_v[15] = imag_x[11] - imag_x[15];
        /* RED */       real_v[3] += real_x[13];
        /* RED */       imag_v[3] += imag_x[13];
        /* RED */       real_v[16] = real_x[12] - real_x[15];
        /* RED */       imag_v[16] = imag_x[12] - imag_x[15];
        /* RED */       real_v[3] += real_x[14];
        /* RED */       imag_v[3] += imag_x[14];
        /* RED */       real_v[17] = real_x[13] - real_x[15];
        /* RED */       imag_v[17] = imag_x[13] - imag_x[15];
        /* RED */       real_v[3] += real_x[15];
        /* RED */       imag_v[3] += imag_x[15];
        /* RED */       real_v[18] = real_x[14] - real_x[15];
        /* RED */       imag_v[18] = imag_x[14] - imag_x[15];
        /* RED */       real_v[4] = real_x[16];
        /* RED */       imag_v[4] = imag_x[16];
        /* RED */       real_v[4] += real_x[17];
        /* RED */       imag_v[4] += imag_x[17];
        /* RED */       real_v[19] = real_x[16] - real_x[20];
        /* RED */       imag_v[19] = imag_x[16] - imag_x[20];
        /* RED */       real_v[4] += real_x[18];
        /* RED */       imag_v[4] += imag_x[18];
        /* RED */       real_v[20] = real_x[17] - real_x[20];
        /* RED */       imag_v[20] = imag_x[17] - imag_x[20];
        /* RED */       real_v[4] += real_x[19];
        /* RED */       imag_v[4] += imag_x[19];
        /* RED */       real_v[21] = real_x[18] - real_x[20];
        /* RED */       imag_v[21] = imag_x[18] - imag_x[20];
        /* RED */       real_v[4] += real_x[20];
        /* RED */       imag_v[4] += imag_x[20];
        /* RED */       real_v[22] = real_x[19] - real_x[20];
        /* RED */       imag_v[22] = imag_x[19] - imag_x[20];
        /* RED */       real_v[5] = real_x[21];
        /* RED */       imag_v[5] = imag_x[21];
        /* RED */       real_v[5] += real_x[22];
        /* RED */       imag_v[5] += imag_x[22];
        /* RED */       real_v[23] = real_x[21] - real_x[25];
        /* RED */       imag_v[23] = imag_x[21] - imag_x[25];
        /* RED */       real_v[5] += real_x[23];
        /* RED */       imag_v[5] += imag_x[23];
        /* RED */       real_v[24] = real_x[22] - real_x[25];
        /* RED */       imag_v[24] = imag_x[22] - imag_x[25];
        /* RED */       real_v[5] += real_x[24];
        /* RED */       imag_v[5] += imag_x[24];
        /* RED */       real_v[25] = real_x[23] - real_x[25];
        /* RED */       imag_v[25] = imag_x[23] - imag_x[25];
        /* RED */       real_v[5] += real_x[25];
        /* RED */       imag_v[5] += imag_x[25];
        /* RED */       real_v[26] = real_x[24] - real_x[25];
        /* RED */       imag_v[26] = imag_x[24] - imag_x[25];
        /* RED */       real_v[6] = real_x[26];
        /* RED */       imag_v[6] = imag_x[26];
        /* RED */       real_v[6] += real_x[27];
        /* RED */       imag_v[6] += imag_x[27];
        /* RED */       real_v[27] = real_x[26] - real_x[30];
        /* RED */       imag_v[27] = imag_x[26] - imag_x[30];
        /* RED */       real_v[6] += real_x[28];
        /* RED */       imag_v[6] += imag_x[28];
        /* RED */       real_v[28] = real_x[27] - real_x[30];
        /* RED */       imag_v[28] = imag_x[27] - imag_x[30];
        /* RED */       real_v[6] += real_x[29];
        /* RED */       imag_v[6] += imag_x[29];
        /* RED */       real_v[29] = real_x[28] - real_x[30];
        /* RED */       imag_v[29] = imag_x[28] - imag_x[30];
        /* RED */       real_v[6] += real_x[30];
        /* RED */       imag_v[6] += imag_x[30];
        /* RED */       real_v[30] = real_x[29] - real_x[30];
        /* RED */       imag_v[30] = imag_x[29] - imag_x[30];
        /* RED */       real_x[1] = real_v[1];
        /* RED */       imag_x[1] = imag_v[1];
        /* RED */       real_x[2] = real_v[2];
        /* RED */       imag_x[2] = imag_v[2];
        /* RED */       real_x[3] = real_v[3];
        /* RED */       imag_x[3] = imag_v[3];
        /* RED */       real_x[4] = real_v[4];
        /* RED */       imag_x[4] = imag_v[4];
        /* RED */       real_x[5] = real_v[5];
        /* RED */       imag_x[5] = imag_v[5];
        /* RED */       real_x[6] = real_v[6];
        /* RED */       imag_x[6] = imag_v[6];
        /* RED */       real_x[7] = real_v[7];
        /* RED */       imag_x[7] = imag_v[7];
        /* RED */       real_x[8] = real_v[8];
        /* RED */       imag_x[8] = imag_v[8];
        /* RED */       real_x[9] = real_v[9];
        /* RED */       imag_x[9] = imag_v[9];
        /* RED */       real_x[10] = real_v[10];
        /* RED */       imag_x[10] = imag_v[10];
        /* RED */       real_x[11] = real_v[11];
        /* RED */       imag_x[11] = imag_v[11];
        /* RED */       real_x[12] = real_v[12];
        /* RED */       imag_x[12] = imag_v[12];
        /* RED */       real_x[13] = real_v[13];
        /* RED */       imag_x[13] = imag_v[13];
        /* RED */       real_x[14] = real_v[14];
        /* RED */       imag_x[14] = imag_v[14];
        /* RED */       real_x[15] = real_v[15];
        /* RED */       imag_x[15] = imag_v[15];
        /* RED */       real_x[16] = real_v[16];
        /* RED */       imag_x[16] = imag_v[16];
        /* RED */       real_x[17] = real_v[17];
        /* RED */       imag_x[17] = imag_v[17];
        /* RED */       real_x[18] = real_v[18];
        /* RED */       imag_x[18] = imag_v[18];
        /* RED */       real_x[19] = real_v[19];
        /* RED */       imag_x[19] = imag_v[19];
        /* RED */       real_x[20] = real_v[20];
        /* RED */       imag_x[20] = imag_v[20];
        /* RED */       real_x[21] = real_v[21];
        /* RED */       imag_x[21] = imag_v[21];
        /* RED */       real_x[22] = real_v[22];
        /* RED */       imag_x[22] = imag_v[22];
        /* RED */       real_x[23] = real_v[23];
        /* RED */       imag_x[23] = imag_v[23];
        /* RED */       real_x[24] = real_v[24];
        /* RED */       imag_x[24] = imag_v[24];
        /* RED */       real_x[25] = real_v[25];
        /* RED */       imag_x[25] = imag_v[25];
        /* RED */       real_x[26] = real_v[26];
        /* RED */       imag_x[26] = imag_v[26];
        /* RED */       real_x[27] = real_v[27];
        /* RED */       imag_x[27] = imag_v[27];
        /* RED */       real_x[28] = real_v[28];
        /* RED */       imag_x[28] = imag_v[28];
        /* RED */       real_x[29] = real_v[29];
        /* RED */       imag_x[29] = imag_v[29];
        /* RED */       real_x[30] = real_v[30];
        /* RED */       imag_x[30] = imag_v[30];
        /* RED */       /*  RED exit */
      /* KRED */       /*  KRED exit  */
        /* ADD  */       real_y[0] = real_x[0] + real_x[1];
        /* ADD  */       imag_y[0] = imag_x[0] + imag_x[1];
        /* MULre  */       real_y[1] = real_x[1] * u[0];
        /* MULre  */       imag_y[1] = imag_x[1] * u[0];
        /* MULim  */       real_y[2] = -1 * imag_x[2] * u[1];
        /* MULim  */       imag_y[2] = real_x[2] * u[1];
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[3];
        /* D2  */       imag_v[0] = imag_x[3];
        /* D2  */       real_v[1] = real_x[4];
        /* D2  */       imag_v[1] = imag_x[4];
        /* D2  */       real_v[2] = real_x[3] + real_x[4];
        /* D2  */       imag_v[2] = imag_x[3] + imag_x[4];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[2];
        /* REAL  */       imag_v[0] = imag_v[0] * u[2];
        /* REAL  */       real_v[1] = real_v[1] * u[3];
        /* REAL  */       imag_v[1] = imag_v[1] * u[3];
        /* REAL  */       real_v[2] = real_v[2] * u[4];
        /* REAL  */       imag_v[2] = imag_v[2] * u[4];
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[3] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[3] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[4] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[4] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 1 */
        /* D2  */       real_v[0] = real_x[5];
        /* D2  */       imag_v[0] = imag_x[5];
        /* D2  */       real_v[1] = real_x[6];
        /* D2  */       imag_v[1] = imag_x[6];
        /* D2  */       real_v[2] = real_x[5] + real_x[6];
        /* D2  */       imag_v[2] = imag_x[5] + imag_x[6];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[5];
        /* IMAG  */       imag_t = real_v[0] * u[5];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[6];
        /* IMAG  */       imag_t = real_v[1] * u[6];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[7];
        /* IMAG  */       imag_t = real_v[2] * u[7];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 1  */
        /* D2t  */       real_y[5] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[5] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[6] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[6] = imag_v[1] + imag_v[2];
        /* ID2It */     /* Exit */
     /* MARKER */
        /* ID2I */     /* Entry m = 1  n = 2 */
        /* D2  */       real_v1[0] = real_x[7];
        /* D2  */       imag_v1[0] = imag_x[7];
        /* D2  */       real_v1[2] = real_x[9];
        /* D2  */       imag_v1[2] = imag_x[9];
        /* D2  */       real_v1[4] = real_x[7] + real_x[9];
        /* D2  */       imag_v1[4] = imag_x[7] + imag_x[9];
        /* D2  */       real_v1[1] = real_x[8];
        /* D2  */       imag_v1[1] = imag_x[8];
        /* D2  */       real_v1[3] = real_x[10];
        /* D2  */       imag_v1[3] = imag_x[10];
        /* D2  */       real_v1[5] = real_x[8] + real_x[10];
        /* D2  */       imag_v1[5] = imag_x[8] + imag_x[10];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 1 */
        /* D2  */       real_v[0] = real_v1[0];
        /* D2  */       imag_v[0] = imag_v1[0];
        /* D2  */       real_v[1] = real_v1[1];
        /* D2  */       imag_v[1] = imag_v1[1];
        /* D2  */       real_v[2] = real_v1[0] + real_v1[1];
        /* D2  */       imag_v[2] = imag_v1[0] + imag_v1[1];
        /* D2  */       real_v[3] = real_v1[2];
        /* D2  */       imag_v[3] = imag_v1[2];
        /* D2  */       real_v[4] = real_v1[3];
        /* D2  */       imag_v[4] = imag_v1[3];
        /* D2  */       real_v[5] = real_v1[2] + real_v1[3];
        /* D2  */       imag_v[5] = imag_v1[2] + imag_v1[3];
        /* D2  */       real_v[6] = real_v1[4];
        /* D2  */       imag_v[6] = imag_v1[4];
        /* D2  */       real_v[7] = real_v1[5];
        /* D2  */       imag_v[7] = imag_v1[5];
        /* D2  */       real_v[8] = real_v1[4] + real_v1[5];
        /* D2  */       imag_v[8] = imag_v1[4] + imag_v1[5];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[8];
        /* REAL  */       imag_v[0] = imag_v[0] * u[8];
        /* REAL  */       real_v[1] = real_v[1] * u[9];
        /* REAL  */       imag_v[1] = imag_v[1] * u[9];
        /* REAL  */       real_v[2] = real_v[2] * u[10];
        /* REAL  */       imag_v[2] = imag_v[2] * u[10];
        /* REAL  */       real_v[3] = real_v[3] * u[11];
        /* REAL  */       imag_v[3] = imag_v[3] * u[11];
        /* REAL  */       real_v[4] = real_v[4] * u[12];
        /* REAL  */       imag_v[4] = imag_v[4] * u[12];
        /* REAL  */       real_v[5] = real_v[5] * u[13];
        /* REAL  */       imag_v[5] = imag_v[5] * u[13];
        /* REAL  */       real_v[6] = real_v[6] * u[14];
        /* REAL  */       imag_v[6] = imag_v[6] * u[14];
        /* REAL  */       real_v[7] = real_v[7] * u[15];
        /* REAL  */       imag_v[7] = imag_v[7] * u[15];
        /* REAL  */       real_v[8] = real_v[8] * u[16];
        /* REAL  */       imag_v[8] = imag_v[8] * u[16];
        /* ID2It */     /* Entry m = 1  n = 3  */
        /* D2t  */       real_v1[0] = real_v[0] + real_v[6];
        /* D2t  */       imag_v1[0] = imag_v[0] + imag_v[6];
        /* D2t  */       real_v1[3] = real_v[3] + real_v[6];
        /* D2t  */       imag_v1[3] = imag_v[3] + imag_v[6];
        /* D2t  */       real_v1[1] = real_v[1] + real_v[7];
        /* D2t  */       imag_v1[1] = imag_v[1] + imag_v[7];
        /* D2t  */       real_v1[4] = real_v[4] + real_v[7];
        /* D2t  */       imag_v1[4] = imag_v[4] + imag_v[7];
        /* D2t  */       real_v1[2] = real_v[2] + real_v[8];
        /* D2t  */       imag_v1[2] = imag_v[2] + imag_v[8];
        /* D2t  */       real_v1[5] = real_v[5] + real_v[8];
        /* D2t  */       imag_v1[5] = imag_v[5] + imag_v[8];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 1  */
        /* D2t  */       real_y[7] = real_v1[0] + real_v1[2];
        /* D2t  */       imag_y[7] = imag_v1[0] + imag_v1[2];
        /* D2t  */       real_y[8] = real_v1[1] + real_v1[2];
        /* D2t  */       imag_y[8] = imag_v1[1] + imag_v1[2];
        /* D2t  */       real_y[9] = real_v1[3] + real_v1[5];
        /* D2t  */       imag_y[9] = imag_v1[3] + imag_v1[5];
        /* D2t  */       real_y[10] = real_v1[4] + real_v1[5];
        /* D2t  */       imag_y[10] = imag_v1[4] + imag_v1[5];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 2 */
        /* D2  */       real_v1[0] = real_x[11];
        /* D2  */       imag_v1[0] = imag_x[11];
        /* D2  */       real_v1[2] = real_x[13];
        /* D2  */       imag_v1[2] = imag_x[13];
        /* D2  */       real_v1[4] = real_x[11] + real_x[13];
        /* D2  */       imag_v1[4] = imag_x[11] + imag_x[13];
        /* D2  */       real_v1[1] = real_x[12];
        /* D2  */       imag_v1[1] = imag_x[12];
        /* D2  */       real_v1[3] = real_x[14];
        /* D2  */       imag_v1[3] = imag_x[14];
        /* D2  */       real_v1[5] = real_x[12] + real_x[14];
        /* D2  */       imag_v1[5] = imag_x[12] + imag_x[14];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 1 */
        /* D2  */       real_v[0] = real_v1[0];
        /* D2  */       imag_v[0] = imag_v1[0];
        /* D2  */       real_v[1] = real_v1[1];
        /* D2  */       imag_v[1] = imag_v1[1];
        /* D2  */       real_v[2] = real_v1[0] + real_v1[1];
        /* D2  */       imag_v[2] = imag_v1[0] + imag_v1[1];
        /* D2  */       real_v[3] = real_v1[2];
        /* D2  */       imag_v[3] = imag_v1[2];
        /* D2  */       real_v[4] = real_v1[3];
        /* D2  */       imag_v[4] = imag_v1[3];
        /* D2  */       real_v[5] = real_v1[2] + real_v1[3];
        /* D2  */       imag_v[5] = imag_v1[2] + imag_v1[3];
        /* D2  */       real_v[6] = real_v1[4];
        /* D2  */       imag_v[6] = imag_v1[4];
        /* D2  */       real_v[7] = real_v1[5];
        /* D2  */       imag_v[7] = imag_v1[5];
        /* D2  */       real_v[8] = real_v1[4] + real_v1[5];
        /* D2  */       imag_v[8] = imag_v1[4] + imag_v1[5];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[17];
        /* IMAG  */       imag_t = real_v[0] * u[17];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[18];
        /* IMAG  */       imag_t = real_v[1] * u[18];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[19];
        /* IMAG  */       imag_t = real_v[2] * u[19];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[3] * u[20];
        /* IMAG  */       imag_t = real_v[3] * u[20];
        /* IMAG  */       real_v[3] = real_t;
        /* IMAG  */       imag_v[3] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[4] * u[21];
        /* IMAG  */       imag_t = real_v[4] * u[21];
        /* IMAG  */       real_v[4] = real_t;
        /* IMAG  */       imag_v[4] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[5] * u[22];
        /* IMAG  */       imag_t = real_v[5] * u[22];
        /* IMAG  */       real_v[5] = real_t;
        /* IMAG  */       imag_v[5] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[6] * u[23];
        /* IMAG  */       imag_t = real_v[6] * u[23];
        /* IMAG  */       real_v[6] = real_t;
        /* IMAG  */       imag_v[6] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[7] * u[24];
        /* IMAG  */       imag_t = real_v[7] * u[24];
        /* IMAG  */       real_v[7] = real_t;
        /* IMAG  */       imag_v[7] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[8] * u[25];
        /* IMAG  */       imag_t = real_v[8] * u[25];
        /* IMAG  */       real_v[8] = real_t;
        /* IMAG  */       imag_v[8] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 3  */
        /* D2t  */       real_v1[0] = real_v[0] + real_v[6];
        /* D2t  */       imag_v1[0] = imag_v[0] + imag_v[6];
        /* D2t  */       real_v1[3] = real_v[3] + real_v[6];
        /* D2t  */       imag_v1[3] = imag_v[3] + imag_v[6];
        /* D2t  */       real_v1[1] = real_v[1] + real_v[7];
        /* D2t  */       imag_v1[1] = imag_v[1] + imag_v[7];
        /* D2t  */       real_v1[4] = real_v[4] + real_v[7];
        /* D2t  */       imag_v1[4] = imag_v[4] + imag_v[7];
        /* D2t  */       real_v1[2] = real_v[2] + real_v[8];
        /* D2t  */       imag_v1[2] = imag_v[2] + imag_v[8];
        /* D2t  */       real_v1[5] = real_v[5] + real_v[8];
        /* D2t  */       imag_v1[5] = imag_v[5] + imag_v[8];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 1  */
        /* D2t  */       real_y[11] = real_v1[0] + real_v1[2];
        /* D2t  */       imag_y[11] = imag_v1[0] + imag_v1[2];
        /* D2t  */       real_y[12] = real_v1[1] + real_v1[2];
        /* D2t  */       imag_y[12] = imag_v1[1] + imag_v1[2];
        /* D2t  */       real_y[13] = real_v1[3] + real_v1[5];
        /* D2t  */       imag_y[13] = imag_v1[3] + imag_v1[5];
        /* D2t  */       real_y[14] = real_v1[4] + real_v1[5];
        /* D2t  */       imag_y[14] = imag_v1[4] + imag_v1[5];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 4 */
        /* D2  */       real_v[0] = real_x[15];
        /* D2  */       imag_v[0] = imag_x[15];
        /* D2  */       real_v[4] = real_x[19];
        /* D2  */       imag_v[4] = imag_x[19];
        /* D2  */       real_v[8] = real_x[15] + real_x[19];
        /* D2  */       imag_v[8] = imag_x[15] + imag_x[19];
        /* D2  */       real_v[1] = real_x[16];
        /* D2  */       imag_v[1] = imag_x[16];
        /* D2  */       real_v[5] = real_x[20];
        /* D2  */       imag_v[5] = imag_x[20];
        /* D2  */       real_v[9] = real_x[16] + real_x[20];
        /* D2  */       imag_v[9] = imag_x[16] + imag_x[20];
        /* D2  */       real_v[2] = real_x[17];
        /* D2  */       imag_v[2] = imag_x[17];
        /* D2  */       real_v[6] = real_x[21];
        /* D2  */       imag_v[6] = imag_x[21];
        /* D2  */       real_v[10] = real_x[17] + real_x[21];
        /* D2  */       imag_v[10] = imag_x[17] + imag_x[21];
        /* D2  */       real_v[3] = real_x[18];
        /* D2  */       imag_v[3] = imag_x[18];
        /* D2  */       real_v[7] = real_x[22];
        /* D2  */       imag_v[7] = imag_x[22];
        /* D2  */       real_v[11] = real_x[18] + real_x[22];
        /* D2  */       imag_v[11] = imag_x[18] + imag_x[22];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 2 */
        /* D2  */       real_v1[0] = real_v[0];
        /* D2  */       imag_v1[0] = imag_v[0];
        /* D2  */       real_v1[2] = real_v[2];
        /* D2  */       imag_v1[2] = imag_v[2];
        /* D2  */       real_v1[4] = real_v[0] + real_v[2];
        /* D2  */       imag_v1[4] = imag_v[0] + imag_v[2];
        /* D2  */       real_v1[1] = real_v[1];
        /* D2  */       imag_v1[1] = imag_v[1];
        /* D2  */       real_v1[3] = real_v[3];
        /* D2  */       imag_v1[3] = imag_v[3];
        /* D2  */       real_v1[5] = real_v[1] + real_v[3];
        /* D2  */       imag_v1[5] = imag_v[1] + imag_v[3];
        /* D2  */       real_v1[6] = real_v[4];
        /* D2  */       imag_v1[6] = imag_v[4];
        /* D2  */       real_v1[8] = real_v[6];
        /* D2  */       imag_v1[8] = imag_v[6];
        /* D2  */       real_v1[10] = real_v[4] + real_v[6];
        /* D2  */       imag_v1[10] = imag_v[4] + imag_v[6];
        /* D2  */       real_v1[7] = real_v[5];
        /* D2  */       imag_v1[7] = imag_v[5];
        /* D2  */       real_v1[9] = real_v[7];
        /* D2  */       imag_v1[9] = imag_v[7];
        /* D2  */       real_v1[11] = real_v[5] + real_v[7];
        /* D2  */       imag_v1[11] = imag_v[5] + imag_v[7];
        /* D2  */       real_v1[12] = real_v[8];
        /* D2  */       imag_v1[12] = imag_v[8];
        /* D2  */       real_v1[14] = real_v[10];
        /* D2  */       imag_v1[14] = imag_v[10];
        /* D2  */       real_v1[16] = real_v[8] + real_v[10];
        /* D2  */       imag_v1[16] = imag_v[8] + imag_v[10];
        /* D2  */       real_v1[13] = real_v[9];
        /* D2  */       imag_v1[13] = imag_v[9];
        /* D2  */       real_v1[15] = real_v[11];
        /* D2  */       imag_v1[15] = imag_v[11];
        /* D2  */       real_v1[17] = real_v[9] + real_v[11];
        /* D2  */       imag_v1[17] = imag_v[9] + imag_v[11];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 9  n = 1 */
        /* D2  */       real_v[0] = real_v1[0];
        /* D2  */       imag_v[0] = imag_v1[0];
        /* D2  */       real_v[1] = real_v1[1];
        /* D2  */       imag_v[1] = imag_v1[1];
        /* D2  */       real_v[2] = real_v1[0] + real_v1[1];
        /* D2  */       imag_v[2] = imag_v1[0] + imag_v1[1];
        /* D2  */       real_v[3] = real_v1[2];
        /* D2  */       imag_v[3] = imag_v1[2];
        /* D2  */       real_v[4] = real_v1[3];
        /* D2  */       imag_v[4] = imag_v1[3];
        /* D2  */       real_v[5] = real_v1[2] + real_v1[3];
        /* D2  */       imag_v[5] = imag_v1[2] + imag_v1[3];
        /* D2  */       real_v[6] = real_v1[4];
        /* D2  */       imag_v[6] = imag_v1[4];
        /* D2  */       real_v[7] = real_v1[5];
        /* D2  */       imag_v[7] = imag_v1[5];
        /* D2  */       real_v[8] = real_v1[4] + real_v1[5];
        /* D2  */       imag_v[8] = imag_v1[4] + imag_v1[5];
        /* D2  */       real_v[9] = real_v1[6];
        /* D2  */       imag_v[9] = imag_v1[6];
        /* D2  */       real_v[10] = real_v1[7];
        /* D2  */       imag_v[10] = imag_v1[7];
        /* D2  */       real_v[11] = real_v1[6] + real_v1[7];
        /* D2  */       imag_v[11] = imag_v1[6] + imag_v1[7];
        /* D2  */       real_v[12] = real_v1[8];
        /* D2  */       imag_v[12] = imag_v1[8];
        /* D2  */       real_v[13] = real_v1[9];
        /* D2  */       imag_v[13] = imag_v1[9];
        /* D2  */       real_v[14] = real_v1[8] + real_v1[9];
        /* D2  */       imag_v[14] = imag_v1[8] + imag_v1[9];
        /* D2  */       real_v[15] = real_v1[10];
        /* D2  */       imag_v[15] = imag_v1[10];
        /* D2  */       real_v[16] = real_v1[11];
        /* D2  */       imag_v[16] = imag_v1[11];
        /* D2  */       real_v[17] = real_v1[10] + real_v1[11];
        /* D2  */       imag_v[17] = imag_v1[10] + imag_v1[11];
        /* D2  */       real_v[18] = real_v1[12];
        /* D2  */       imag_v[18] = imag_v1[12];
        /* D2  */       real_v[19] = real_v1[13];
        /* D2  */       imag_v[19] = imag_v1[13];
        /* D2  */       real_v[20] = real_v1[12] + real_v1[13];
        /* D2  */       imag_v[20] = imag_v1[12] + imag_v1[13];
        /* D2  */       real_v[21] = real_v1[14];
        /* D2  */       imag_v[21] = imag_v1[14];
        /* D2  */       real_v[22] = real_v1[15];
        /* D2  */       imag_v[22] = imag_v1[15];
        /* D2  */       real_v[23] = real_v1[14] + real_v1[15];
        /* D2  */       imag_v[23] = imag_v1[14] + imag_v1[15];
        /* D2  */       real_v[24] = real_v1[16];
        /* D2  */       imag_v[24] = imag_v1[16];
        /* D2  */       real_v[25] = real_v1[17];
        /* D2  */       imag_v[25] = imag_v1[17];
        /* D2  */       real_v[26] = real_v1[16] + real_v1[17];
        /* D2  */       imag_v[26] = imag_v1[16] + imag_v1[17];
        /* ID2I */     /* Exit */
        /* REAL  */       real_v[0] = real_v[0] * u[26];
        /* REAL  */       imag_v[0] = imag_v[0] * u[26];
        /* REAL  */       real_v[1] = real_v[1] * u[27];
        /* REAL  */       imag_v[1] = imag_v[1] * u[27];
        /* REAL  */       real_v[2] = real_v[2] * u[28];
        /* REAL  */       imag_v[2] = imag_v[2] * u[28];
        /* REAL  */       real_v[3] = real_v[3] * u[29];
        /* REAL  */       imag_v[3] = imag_v[3] * u[29];
        /* REAL  */       real_v[4] = real_v[4] * u[30];
        /* REAL  */       imag_v[4] = imag_v[4] * u[30];
        /* REAL  */       real_v[5] = real_v[5] * u[31];
        /* REAL  */       imag_v[5] = imag_v[5] * u[31];
        /* REAL  */       real_v[6] = real_v[6] * u[32];
        /* REAL  */       imag_v[6] = imag_v[6] * u[32];
        /* REAL  */       real_v[7] = real_v[7] * u[33];
        /* REAL  */       imag_v[7] = imag_v[7] * u[33];
        /* REAL  */       real_v[8] = real_v[8] * u[34];
        /* REAL  */       imag_v[8] = imag_v[8] * u[34];
        /* REAL  */       real_v[9] = real_v[9] * u[35];
        /* REAL  */       imag_v[9] = imag_v[9] * u[35];
        /* REAL  */       real_v[10] = real_v[10] * u[36];
        /* REAL  */       imag_v[10] = imag_v[10] * u[36];
        /* REAL  */       real_v[11] = real_v[11] * u[37];
        /* REAL  */       imag_v[11] = imag_v[11] * u[37];
        /* REAL  */       real_v[12] = real_v[12] * u[38];
        /* REAL  */       imag_v[12] = imag_v[12] * u[38];
        /* REAL  */       real_v[13] = real_v[13] * u[39];
        /* REAL  */       imag_v[13] = imag_v[13] * u[39];
        /* REAL  */       real_v[14] = real_v[14] * u[40];
        /* REAL  */       imag_v[14] = imag_v[14] * u[40];
        /* REAL  */       real_v[15] = real_v[15] * u[41];
        /* REAL  */       imag_v[15] = imag_v[15] * u[41];
        /* REAL  */       real_v[16] = real_v[16] * u[42];
        /* REAL  */       imag_v[16] = imag_v[16] * u[42];
        /* REAL  */       real_v[17] = real_v[17] * u[43];
        /* REAL  */       imag_v[17] = imag_v[17] * u[43];
        /* REAL  */       real_v[18] = real_v[18] * u[44];
        /* REAL  */       imag_v[18] = imag_v[18] * u[44];
        /* REAL  */       real_v[19] = real_v[19] * u[45];
        /* REAL  */       imag_v[19] = imag_v[19] * u[45];
        /* REAL  */       real_v[20] = real_v[20] * u[46];
        /* REAL  */       imag_v[20] = imag_v[20] * u[46];
        /* REAL  */       real_v[21] = real_v[21] * u[47];
        /* REAL  */       imag_v[21] = imag_v[21] * u[47];
        /* REAL  */       real_v[22] = real_v[22] * u[48];
        /* REAL  */       imag_v[22] = imag_v[22] * u[48];
        /* REAL  */       real_v[23] = real_v[23] * u[49];
        /* REAL  */       imag_v[23] = imag_v[23] * u[49];
        /* REAL  */       real_v[24] = real_v[24] * u[50];
        /* REAL  */       imag_v[24] = imag_v[24] * u[50];
        /* REAL  */       real_v[25] = real_v[25] * u[51];
        /* REAL  */       imag_v[25] = imag_v[25] * u[51];
        /* REAL  */       real_v[26] = real_v[26] * u[52];
        /* REAL  */       imag_v[26] = imag_v[26] * u[52];
        /* ID2It */     /* Entry m = 1  n = 9  */
        /* D2t  */       real_v1[0] = real_v[0] + real_v[18];
        /* D2t  */       imag_v1[0] = imag_v[0] + imag_v[18];
        /* D2t  */       real_v1[9] = real_v[9] + real_v[18];
        /* D2t  */       imag_v1[9] = imag_v[9] + imag_v[18];
        /* D2t  */       real_v1[1] = real_v[1] + real_v[19];
        /* D2t  */       imag_v1[1] = imag_v[1] + imag_v[19];
        /* D2t  */       real_v1[10] = real_v[10] + real_v[19];
        /* D2t  */       imag_v1[10] = imag_v[10] + imag_v[19];
        /* D2t  */       real_v1[2] = real_v[2] + real_v[20];
        /* D2t  */       imag_v1[2] = imag_v[2] + imag_v[20];
        /* D2t  */       real_v1[11] = real_v[11] + real_v[20];
        /* D2t  */       imag_v1[11] = imag_v[11] + imag_v[20];
        /* D2t  */       real_v1[3] = real_v[3] + real_v[21];
        /* D2t  */       imag_v1[3] = imag_v[3] + imag_v[21];
        /* D2t  */       real_v1[12] = real_v[12] + real_v[21];
        /* D2t  */       imag_v1[12] = imag_v[12] + imag_v[21];
        /* D2t  */       real_v1[4] = real_v[4] + real_v[22];
        /* D2t  */       imag_v1[4] = imag_v[4] + imag_v[22];
        /* D2t  */       real_v1[13] = real_v[13] + real_v[22];
        /* D2t  */       imag_v1[13] = imag_v[13] + imag_v[22];
        /* D2t  */       real_v1[5] = real_v[5] + real_v[23];
        /* D2t  */       imag_v1[5] = imag_v[5] + imag_v[23];
        /* D2t  */       real_v1[14] = real_v[14] + real_v[23];
        /* D2t  */       imag_v1[14] = imag_v[14] + imag_v[23];
        /* D2t  */       real_v1[6] = real_v[6] + real_v[24];
        /* D2t  */       imag_v1[6] = imag_v[6] + imag_v[24];
        /* D2t  */       real_v1[15] = real_v[15] + real_v[24];
        /* D2t  */       imag_v1[15] = imag_v[15] + imag_v[24];
        /* D2t  */       real_v1[7] = real_v[7] + real_v[25];
        /* D2t  */       imag_v1[7] = imag_v[7] + imag_v[25];
        /* D2t  */       real_v1[16] = real_v[16] + real_v[25];
        /* D2t  */       imag_v1[16] = imag_v[16] + imag_v[25];
        /* D2t  */       real_v1[8] = real_v[8] + real_v[26];
        /* D2t  */       imag_v1[8] = imag_v[8] + imag_v[26];
        /* D2t  */       real_v1[17] = real_v[17] + real_v[26];
        /* D2t  */       imag_v1[17] = imag_v[17] + imag_v[26];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 3  */
        /* D2t  */       real_v[0] = real_v1[0] + real_v1[6];
        /* D2t  */       imag_v[0] = imag_v1[0] + imag_v1[6];
        /* D2t  */       real_v[3] = real_v1[3] + real_v1[6];
        /* D2t  */       imag_v[3] = imag_v1[3] + imag_v1[6];
        /* D2t  */       real_v[1] = real_v1[1] + real_v1[7];
        /* D2t  */       imag_v[1] = imag_v1[1] + imag_v1[7];
        /* D2t  */       real_v[4] = real_v1[4] + real_v1[7];
        /* D2t  */       imag_v[4] = imag_v1[4] + imag_v1[7];
        /* D2t  */       real_v[2] = real_v1[2] + real_v1[8];
        /* D2t  */       imag_v[2] = imag_v1[2] + imag_v1[8];
        /* D2t  */       real_v[5] = real_v1[5] + real_v1[8];
        /* D2t  */       imag_v[5] = imag_v1[5] + imag_v1[8];
        /* D2t  */       real_v[6] = real_v1[9] + real_v1[15];
        /* D2t  */       imag_v[6] = imag_v1[9] + imag_v1[15];
        /* D2t  */       real_v[9] = real_v1[12] + real_v1[15];
        /* D2t  */       imag_v[9] = imag_v1[12] + imag_v1[15];
        /* D2t  */       real_v[7] = real_v1[10] + real_v1[16];
        /* D2t  */       imag_v[7] = imag_v1[10] + imag_v1[16];
        /* D2t  */       real_v[10] = real_v1[13] + real_v1[16];
        /* D2t  */       imag_v[10] = imag_v1[13] + imag_v1[16];
        /* D2t  */       real_v[8] = real_v1[11] + real_v1[17];
        /* D2t  */       imag_v[8] = imag_v1[11] + imag_v1[17];
        /* D2t  */       real_v[11] = real_v1[14] + real_v1[17];
        /* D2t  */       imag_v[11] = imag_v1[14] + imag_v1[17];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 4  n = 1  */
        /* D2t  */       real_y[15] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[15] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[16] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[16] = imag_v[1] + imag_v[2];
        /* D2t  */       real_y[17] = real_v[3] + real_v[5];
        /* D2t  */       imag_y[17] = imag_v[3] + imag_v[5];
        /* D2t  */       real_y[18] = real_v[4] + real_v[5];
        /* D2t  */       imag_y[18] = imag_v[4] + imag_v[5];
        /* D2t  */       real_y[19] = real_v[6] + real_v[8];
        /* D2t  */       imag_y[19] = imag_v[6] + imag_v[8];
        /* D2t  */       real_y[20] = real_v[7] + real_v[8];
        /* D2t  */       imag_y[20] = imag_v[7] + imag_v[8];
        /* D2t  */       real_y[21] = real_v[9] + real_v[11];
        /* D2t  */       imag_y[21] = imag_v[9] + imag_v[11];
        /* D2t  */       real_y[22] = real_v[10] + real_v[11];
        /* D2t  */       imag_y[22] = imag_v[10] + imag_v[11];
        /* ID2It */     /* Exit */
        /* ID2I */     /* Entry m = 1  n = 4 */
        /* D2  */       real_v[0] = real_x[23];
        /* D2  */       imag_v[0] = imag_x[23];
        /* D2  */       real_v[4] = real_x[27];
        /* D2  */       imag_v[4] = imag_x[27];
        /* D2  */       real_v[8] = real_x[23] + real_x[27];
        /* D2  */       imag_v[8] = imag_x[23] + imag_x[27];
        /* D2  */       real_v[1] = real_x[24];
        /* D2  */       imag_v[1] = imag_x[24];
        /* D2  */       real_v[5] = real_x[28];
        /* D2  */       imag_v[5] = imag_x[28];
        /* D2  */       real_v[9] = real_x[24] + real_x[28];
        /* D2  */       imag_v[9] = imag_x[24] + imag_x[28];
        /* D2  */       real_v[2] = real_x[25];
        /* D2  */       imag_v[2] = imag_x[25];
        /* D2  */       real_v[6] = real_x[29];
        /* D2  */       imag_v[6] = imag_x[29];
        /* D2  */       real_v[10] = real_x[25] + real_x[29];
        /* D2  */       imag_v[10] = imag_x[25] + imag_x[29];
        /* D2  */       real_v[3] = real_x[26];
        /* D2  */       imag_v[3] = imag_x[26];
        /* D2  */       real_v[7] = real_x[30];
        /* D2  */       imag_v[7] = imag_x[30];
        /* D2  */       real_v[11] = real_x[26] + real_x[30];
        /* D2  */       imag_v[11] = imag_x[26] + imag_x[30];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 3  n = 2 */
        /* D2  */       real_v1[0] = real_v[0];
        /* D2  */       imag_v1[0] = imag_v[0];
        /* D2  */       real_v1[2] = real_v[2];
        /* D2  */       imag_v1[2] = imag_v[2];
        /* D2  */       real_v1[4] = real_v[0] + real_v[2];
        /* D2  */       imag_v1[4] = imag_v[0] + imag_v[2];
        /* D2  */       real_v1[1] = real_v[1];
        /* D2  */       imag_v1[1] = imag_v[1];
        /* D2  */       real_v1[3] = real_v[3];
        /* D2  */       imag_v1[3] = imag_v[3];
        /* D2  */       real_v1[5] = real_v[1] + real_v[3];
        /* D2  */       imag_v1[5] = imag_v[1] + imag_v[3];
        /* D2  */       real_v1[6] = real_v[4];
        /* D2  */       imag_v1[6] = imag_v[4];
        /* D2  */       real_v1[8] = real_v[6];
        /* D2  */       imag_v1[8] = imag_v[6];
        /* D2  */       real_v1[10] = real_v[4] + real_v[6];
        /* D2  */       imag_v1[10] = imag_v[4] + imag_v[6];
        /* D2  */       real_v1[7] = real_v[5];
        /* D2  */       imag_v1[7] = imag_v[5];
        /* D2  */       real_v1[9] = real_v[7];
        /* D2  */       imag_v1[9] = imag_v[7];
        /* D2  */       real_v1[11] = real_v[5] + real_v[7];
        /* D2  */       imag_v1[11] = imag_v[5] + imag_v[7];
        /* D2  */       real_v1[12] = real_v[8];
        /* D2  */       imag_v1[12] = imag_v[8];
        /* D2  */       real_v1[14] = real_v[10];
        /* D2  */       imag_v1[14] = imag_v[10];
        /* D2  */       real_v1[16] = real_v[8] + real_v[10];
        /* D2  */       imag_v1[16] = imag_v[8] + imag_v[10];
        /* D2  */       real_v1[13] = real_v[9];
        /* D2  */       imag_v1[13] = imag_v[9];
        /* D2  */       real_v1[15] = real_v[11];
        /* D2  */       imag_v1[15] = imag_v[11];
        /* D2  */       real_v1[17] = real_v[9] + real_v[11];
        /* D2  */       imag_v1[17] = imag_v[9] + imag_v[11];
        /* ID2I */     /* Exit */
        /* ID2I */     /* Entry m = 9  n = 1 */
        /* D2  */       real_v[0] = real_v1[0];
        /* D2  */       imag_v[0] = imag_v1[0];
        /* D2  */       real_v[1] = real_v1[1];
        /* D2  */       imag_v[1] = imag_v1[1];
        /* D2  */       real_v[2] = real_v1[0] + real_v1[1];
        /* D2  */       imag_v[2] = imag_v1[0] + imag_v1[1];
        /* D2  */       real_v[3] = real_v1[2];
        /* D2  */       imag_v[3] = imag_v1[2];
        /* D2  */       real_v[4] = real_v1[3];
        /* D2  */       imag_v[4] = imag_v1[3];
        /* D2  */       real_v[5] = real_v1[2] + real_v1[3];
        /* D2  */       imag_v[5] = imag_v1[2] + imag_v1[3];
        /* D2  */       real_v[6] = real_v1[4];
        /* D2  */       imag_v[6] = imag_v1[4];
        /* D2  */       real_v[7] = real_v1[5];
        /* D2  */       imag_v[7] = imag_v1[5];
        /* D2  */       real_v[8] = real_v1[4] + real_v1[5];
        /* D2  */       imag_v[8] = imag_v1[4] + imag_v1[5];
        /* D2  */       real_v[9] = real_v1[6];
        /* D2  */       imag_v[9] = imag_v1[6];
        /* D2  */       real_v[10] = real_v1[7];
        /* D2  */       imag_v[10] = imag_v1[7];
        /* D2  */       real_v[11] = real_v1[6] + real_v1[7];
        /* D2  */       imag_v[11] = imag_v1[6] + imag_v1[7];
        /* D2  */       real_v[12] = real_v1[8];
        /* D2  */       imag_v[12] = imag_v1[8];
        /* D2  */       real_v[13] = real_v1[9];
        /* D2  */       imag_v[13] = imag_v1[9];
        /* D2  */       real_v[14] = real_v1[8] + real_v1[9];
        /* D2  */       imag_v[14] = imag_v1[8] + imag_v1[9];
        /* D2  */       real_v[15] = real_v1[10];
        /* D2  */       imag_v[15] = imag_v1[10];
        /* D2  */       real_v[16] = real_v1[11];
        /* D2  */       imag_v[16] = imag_v1[11];
        /* D2  */       real_v[17] = real_v1[10] + real_v1[11];
        /* D2  */       imag_v[17] = imag_v1[10] + imag_v1[11];
        /* D2  */       real_v[18] = real_v1[12];
        /* D2  */       imag_v[18] = imag_v1[12];
        /* D2  */       real_v[19] = real_v1[13];
        /* D2  */       imag_v[19] = imag_v1[13];
        /* D2  */       real_v[20] = real_v1[12] + real_v1[13];
        /* D2  */       imag_v[20] = imag_v1[12] + imag_v1[13];
        /* D2  */       real_v[21] = real_v1[14];
        /* D2  */       imag_v[21] = imag_v1[14];
        /* D2  */       real_v[22] = real_v1[15];
        /* D2  */       imag_v[22] = imag_v1[15];
        /* D2  */       real_v[23] = real_v1[14] + real_v1[15];
        /* D2  */       imag_v[23] = imag_v1[14] + imag_v1[15];
        /* D2  */       real_v[24] = real_v1[16];
        /* D2  */       imag_v[24] = imag_v1[16];
        /* D2  */       real_v[25] = real_v1[17];
        /* D2  */       imag_v[25] = imag_v1[17];
        /* D2  */       real_v[26] = real_v1[16] + real_v1[17];
        /* D2  */       imag_v[26] = imag_v1[16] + imag_v1[17];
        /* ID2I */     /* Exit */
        /* IMAG  */       real_t = -1.0 * imag_v[0] * u[53];
        /* IMAG  */       imag_t = real_v[0] * u[53];
        /* IMAG  */       real_v[0] = real_t;
        /* IMAG  */       imag_v[0] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[1] * u[54];
        /* IMAG  */       imag_t = real_v[1] * u[54];
        /* IMAG  */       real_v[1] = real_t;
        /* IMAG  */       imag_v[1] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[2] * u[55];
        /* IMAG  */       imag_t = real_v[2] * u[55];
        /* IMAG  */       real_v[2] = real_t;
        /* IMAG  */       imag_v[2] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[3] * u[56];
        /* IMAG  */       imag_t = real_v[3] * u[56];
        /* IMAG  */       real_v[3] = real_t;
        /* IMAG  */       imag_v[3] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[4] * u[57];
        /* IMAG  */       imag_t = real_v[4] * u[57];
        /* IMAG  */       real_v[4] = real_t;
        /* IMAG  */       imag_v[4] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[5] * u[58];
        /* IMAG  */       imag_t = real_v[5] * u[58];
        /* IMAG  */       real_v[5] = real_t;
        /* IMAG  */       imag_v[5] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[6] * u[59];
        /* IMAG  */       imag_t = real_v[6] * u[59];
        /* IMAG  */       real_v[6] = real_t;
        /* IMAG  */       imag_v[6] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[7] * u[60];
        /* IMAG  */       imag_t = real_v[7] * u[60];
        /* IMAG  */       real_v[7] = real_t;
        /* IMAG  */       imag_v[7] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[8] * u[61];
        /* IMAG  */       imag_t = real_v[8] * u[61];
        /* IMAG  */       real_v[8] = real_t;
        /* IMAG  */       imag_v[8] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[9] * u[62];
        /* IMAG  */       imag_t = real_v[9] * u[62];
        /* IMAG  */       real_v[9] = real_t;
        /* IMAG  */       imag_v[9] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[10] * u[63];
        /* IMAG  */       imag_t = real_v[10] * u[63];
        /* IMAG  */       real_v[10] = real_t;
        /* IMAG  */       imag_v[10] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[11] * u[64];
        /* IMAG  */       imag_t = real_v[11] * u[64];
        /* IMAG  */       real_v[11] = real_t;
        /* IMAG  */       imag_v[11] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[12] * u[65];
        /* IMAG  */       imag_t = real_v[12] * u[65];
        /* IMAG  */       real_v[12] = real_t;
        /* IMAG  */       imag_v[12] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[13] * u[66];
        /* IMAG  */       imag_t = real_v[13] * u[66];
        /* IMAG  */       real_v[13] = real_t;
        /* IMAG  */       imag_v[13] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[14] * u[67];
        /* IMAG  */       imag_t = real_v[14] * u[67];
        /* IMAG  */       real_v[14] = real_t;
        /* IMAG  */       imag_v[14] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[15] * u[68];
        /* IMAG  */       imag_t = real_v[15] * u[68];
        /* IMAG  */       real_v[15] = real_t;
        /* IMAG  */       imag_v[15] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[16] * u[69];
        /* IMAG  */       imag_t = real_v[16] * u[69];
        /* IMAG  */       real_v[16] = real_t;
        /* IMAG  */       imag_v[16] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[17] * u[70];
        /* IMAG  */       imag_t = real_v[17] * u[70];
        /* IMAG  */       real_v[17] = real_t;
        /* IMAG  */       imag_v[17] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[18] * u[71];
        /* IMAG  */       imag_t = real_v[18] * u[71];
        /* IMAG  */       real_v[18] = real_t;
        /* IMAG  */       imag_v[18] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[19] * u[72];
        /* IMAG  */       imag_t = real_v[19] * u[72];
        /* IMAG  */       real_v[19] = real_t;
        /* IMAG  */       imag_v[19] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[20] * u[73];
        /* IMAG  */       imag_t = real_v[20] * u[73];
        /* IMAG  */       real_v[20] = real_t;
        /* IMAG  */       imag_v[20] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[21] * u[74];
        /* IMAG  */       imag_t = real_v[21] * u[74];
        /* IMAG  */       real_v[21] = real_t;
        /* IMAG  */       imag_v[21] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[22] * u[75];
        /* IMAG  */       imag_t = real_v[22] * u[75];
        /* IMAG  */       real_v[22] = real_t;
        /* IMAG  */       imag_v[22] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[23] * u[76];
        /* IMAG  */       imag_t = real_v[23] * u[76];
        /* IMAG  */       real_v[23] = real_t;
        /* IMAG  */       imag_v[23] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[24] * u[77];
        /* IMAG  */       imag_t = real_v[24] * u[77];
        /* IMAG  */       real_v[24] = real_t;
        /* IMAG  */       imag_v[24] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[25] * u[78];
        /* IMAG  */       imag_t = real_v[25] * u[78];
        /* IMAG  */       real_v[25] = real_t;
        /* IMAG  */       imag_v[25] = imag_t;
        /* IMAG  */       real_t = -1.0 * imag_v[26] * u[79];
        /* IMAG  */       imag_t = real_v[26] * u[79];
        /* IMAG  */       real_v[26] = real_t;
        /* IMAG  */       imag_v[26] = imag_t;
        /* ID2It */     /* Entry m = 1  n = 9  */
        /* D2t  */       real_v1[0] = real_v[0] + real_v[18];
        /* D2t  */       imag_v1[0] = imag_v[0] + imag_v[18];
        /* D2t  */       real_v1[9] = real_v[9] + real_v[18];
        /* D2t  */       imag_v1[9] = imag_v[9] + imag_v[18];
        /* D2t  */       real_v1[1] = real_v[1] + real_v[19];
        /* D2t  */       imag_v1[1] = imag_v[1] + imag_v[19];
        /* D2t  */       real_v1[10] = real_v[10] + real_v[19];
        /* D2t  */       imag_v1[10] = imag_v[10] + imag_v[19];
        /* D2t  */       real_v1[2] = real_v[2] + real_v[20];
        /* D2t  */       imag_v1[2] = imag_v[2] + imag_v[20];
        /* D2t  */       real_v1[11] = real_v[11] + real_v[20];
        /* D2t  */       imag_v1[11] = imag_v[11] + imag_v[20];
        /* D2t  */       real_v1[3] = real_v[3] + real_v[21];
        /* D2t  */       imag_v1[3] = imag_v[3] + imag_v[21];
        /* D2t  */       real_v1[12] = real_v[12] + real_v[21];
        /* D2t  */       imag_v1[12] = imag_v[12] + imag_v[21];
        /* D2t  */       real_v1[4] = real_v[4] + real_v[22];
        /* D2t  */       imag_v1[4] = imag_v[4] + imag_v[22];
        /* D2t  */       real_v1[13] = real_v[13] + real_v[22];
        /* D2t  */       imag_v1[13] = imag_v[13] + imag_v[22];
        /* D2t  */       real_v1[5] = real_v[5] + real_v[23];
        /* D2t  */       imag_v1[5] = imag_v[5] + imag_v[23];
        /* D2t  */       real_v1[14] = real_v[14] + real_v[23];
        /* D2t  */       imag_v1[14] = imag_v[14] + imag_v[23];
        /* D2t  */       real_v1[6] = real_v[6] + real_v[24];
        /* D2t  */       imag_v1[6] = imag_v[6] + imag_v[24];
        /* D2t  */       real_v1[15] = real_v[15] + real_v[24];
        /* D2t  */       imag_v1[15] = imag_v[15] + imag_v[24];
        /* D2t  */       real_v1[7] = real_v[7] + real_v[25];
        /* D2t  */       imag_v1[7] = imag_v[7] + imag_v[25];
        /* D2t  */       real_v1[16] = real_v[16] + real_v[25];
        /* D2t  */       imag_v1[16] = imag_v[16] + imag_v[25];
        /* D2t  */       real_v1[8] = real_v[8] + real_v[26];
        /* D2t  */       imag_v1[8] = imag_v[8] + imag_v[26];
        /* D2t  */       real_v1[17] = real_v[17] + real_v[26];
        /* D2t  */       imag_v1[17] = imag_v[17] + imag_v[26];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 2  n = 3  */
        /* D2t  */       real_v[0] = real_v1[0] + real_v1[6];
        /* D2t  */       imag_v[0] = imag_v1[0] + imag_v1[6];
        /* D2t  */       real_v[3] = real_v1[3] + real_v1[6];
        /* D2t  */       imag_v[3] = imag_v1[3] + imag_v1[6];
        /* D2t  */       real_v[1] = real_v1[1] + real_v1[7];
        /* D2t  */       imag_v[1] = imag_v1[1] + imag_v1[7];
        /* D2t  */       real_v[4] = real_v1[4] + real_v1[7];
        /* D2t  */       imag_v[4] = imag_v1[4] + imag_v1[7];
        /* D2t  */       real_v[2] = real_v1[2] + real_v1[8];
        /* D2t  */       imag_v[2] = imag_v1[2] + imag_v1[8];
        /* D2t  */       real_v[5] = real_v1[5] + real_v1[8];
        /* D2t  */       imag_v[5] = imag_v1[5] + imag_v1[8];
        /* D2t  */       real_v[6] = real_v1[9] + real_v1[15];
        /* D2t  */       imag_v[6] = imag_v1[9] + imag_v1[15];
        /* D2t  */       real_v[9] = real_v1[12] + real_v1[15];
        /* D2t  */       imag_v[9] = imag_v1[12] + imag_v1[15];
        /* D2t  */       real_v[7] = real_v1[10] + real_v1[16];
        /* D2t  */       imag_v[7] = imag_v1[10] + imag_v1[16];
        /* D2t  */       real_v[10] = real_v1[13] + real_v1[16];
        /* D2t  */       imag_v[10] = imag_v1[13] + imag_v1[16];
        /* D2t  */       real_v[8] = real_v1[11] + real_v1[17];
        /* D2t  */       imag_v[8] = imag_v1[11] + imag_v1[17];
        /* D2t  */       real_v[11] = real_v1[14] + real_v1[17];
        /* D2t  */       imag_v[11] = imag_v1[14] + imag_v1[17];
        /* ID2It */     /* Exit */
        /* ID2It */     /* Entry m = 4  n = 1  */
        /* D2t  */       real_y[23] = real_v[0] + real_v[2];
        /* D2t  */       imag_y[23] = imag_v[0] + imag_v[2];
        /* D2t  */       real_y[24] = real_v[1] + real_v[2];
        /* D2t  */       imag_y[24] = imag_v[1] + imag_v[2];
        /* D2t  */       real_y[25] = real_v[3] + real_v[5];
        /* D2t  */       imag_y[25] = imag_v[3] + imag_v[5];
        /* D2t  */       real_y[26] = real_v[4] + real_v[5];
        /* D2t  */       imag_y[26] = imag_v[4] + imag_v[5];
        /* D2t  */       real_y[27] = real_v[6] + real_v[8];
        /* D2t  */       imag_y[27] = imag_v[6] + imag_v[8];
        /* D2t  */       real_y[28] = real_v[7] + real_v[8];
        /* D2t  */       imag_y[28] = imag_v[7] + imag_v[8];
        /* D2t  */       real_y[29] = real_v[9] + real_v[11];
        /* D2t  */       imag_y[29] = imag_v[9] + imag_v[11];
        /* D2t  */       real_y[30] = real_v[10] + real_v[11];
        /* D2t  */       imag_y[30] = imag_v[10] + imag_v[11];
        /* ID2It */     /* Exit */
        /* ADD  */       real_y[1] = real_y[0] + real_y[1];
        /* ADD  */       imag_y[1] = imag_y[0] + imag_y[1];
        /* tKRED */       /* tKRED entry */
          /* tRED */       /*  tRED entry p 5 a 6 c 1  */
        /* tRED */       real_v[5] = real_y[1];
        /* tRED */       imag_v[5] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[7];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[7];
        /* tRED */       real_v[5] -= real_y[7];
        /* tRED */       imag_v[5] -= imag_y[7];
        /* tRED */       real_v[2] = real_y[1] + real_y[8];
        /* tRED */       imag_v[2] = imag_y[1] + imag_y[8];
        /* tRED */       real_v[5] -= real_y[8];
        /* tRED */       imag_v[5] -= imag_y[8];
        /* tRED */       real_v[3] = real_y[1] + real_y[9];
        /* tRED */       imag_v[3] = imag_y[1] + imag_y[9];
        /* tRED */       real_v[5] -= real_y[9];
        /* tRED */       imag_v[5] -= imag_y[9];
        /* tRED */       real_v[4] = real_y[1] + real_y[10];
        /* tRED */       imag_v[4] = imag_y[1] + imag_y[10];
        /* tRED */       real_v[5] -= real_y[10];
        /* tRED */       imag_v[5] -= imag_y[10];
        /* tRED */       real_v[10] = real_y[2];
        /* tRED */       imag_v[10] = imag_y[2];
        /* tRED */       real_v[6] = real_y[2] + real_y[11];
        /* tRED */       imag_v[6] = imag_y[2] + imag_y[11];
        /* tRED */       real_v[10] -= real_y[11];
        /* tRED */       imag_v[10] -= imag_y[11];
        /* tRED */       real_v[7] = real_y[2] + real_y[12];
        /* tRED */       imag_v[7] = imag_y[2] + imag_y[12];
        /* tRED */       real_v[10] -= real_y[12];
        /* tRED */       imag_v[10] -= imag_y[12];
        /* tRED */       real_v[8] = real_y[2] + real_y[13];
        /* tRED */       imag_v[8] = imag_y[2] + imag_y[13];
        /* tRED */       real_v[10] -= real_y[13];
        /* tRED */       imag_v[10] -= imag_y[13];
        /* tRED */       real_v[9] = real_y[2] + real_y[14];
        /* tRED */       imag_v[9] = imag_y[2] + imag_y[14];
        /* tRED */       real_v[10] -= real_y[14];
        /* tRED */       imag_v[10] -= imag_y[14];
        /* tRED */       real_v[15] = real_y[3];
        /* tRED */       imag_v[15] = imag_y[3];
        /* tRED */       real_v[11] = real_y[3] + real_y[15];
        /* tRED */       imag_v[11] = imag_y[3] + imag_y[15];
        /* tRED */       real_v[15] -= real_y[15];
        /* tRED */       imag_v[15] -= imag_y[15];
        /* tRED */       real_v[12] = real_y[3] + real_y[16];
        /* tRED */       imag_v[12] = imag_y[3] + imag_y[16];
        /* tRED */       real_v[15] -= real_y[16];
        /* tRED */       imag_v[15] -= imag_y[16];
        /* tRED */       real_v[13] = real_y[3] + real_y[17];
        /* tRED */       imag_v[13] = imag_y[3] + imag_y[17];
        /* tRED */       real_v[15] -= real_y[17];
        /* tRED */       imag_v[15] -= imag_y[17];
        /* tRED */       real_v[14] = real_y[3] + real_y[18];
        /* tRED */       imag_v[14] = imag_y[3] + imag_y[18];
        /* tRED */       real_v[15] -= real_y[18];
        /* tRED */       imag_v[15] -= imag_y[18];
        /* tRED */       real_v[20] = real_y[4];
        /* tRED */       imag_v[20] = imag_y[4];
        /* tRED */       real_v[16] = real_y[4] + real_y[19];
        /* tRED */       imag_v[16] = imag_y[4] + imag_y[19];
        /* tRED */       real_v[20] -= real_y[19];
        /* tRED */       imag_v[20] -= imag_y[19];
        /* tRED */       real_v[17] = real_y[4] + real_y[20];
        /* tRED */       imag_v[17] = imag_y[4] + imag_y[20];
        /* tRED */       real_v[20] -= real_y[20];
        /* tRED */       imag_v[20] -= imag_y[20];
        /* tRED */       real_v[18] = real_y[4] + real_y[21];
        /* tRED */       imag_v[18] = imag_y[4] + imag_y[21];
        /* tRED */       real_v[20] -= real_y[21];
        /* tRED */       imag_v[20] -= imag_y[21];
        /* tRED */       real_v[19] = real_y[4] + real_y[22];
        /* tRED */       imag_v[19] = imag_y[4] + imag_y[22];
        /* tRED */       real_v[20] -= real_y[22];
        /* tRED */       imag_v[20] -= imag_y[22];
        /* tRED */       real_v[25] = real_y[5];
        /* tRED */       imag_v[25] = imag_y[5];
        /* tRED */       real_v[21] = real_y[5] + real_y[23];
        /* tRED */       imag_v[21] = imag_y[5] + imag_y[23];
        /* tRED */       real_v[25] -= real_y[23];
        /* tRED */       imag_v[25] -= imag_y[23];
        /* tRED */       real_v[22] = real_y[5] + real_y[24];
        /* tRED */       imag_v[22] = imag_y[5] + imag_y[24];
        /* tRED */       real_v[25] -= real_y[24];
        /* tRED */       imag_v[25] -= imag_y[24];
        /* tRED */       real_v[23] = real_y[5] + real_y[25];
        /* tRED */       imag_v[23] = imag_y[5] + imag_y[25];
        /* tRED */       real_v[25] -= real_y[25];
        /* tRED */       imag_v[25] -= imag_y[25];
        /* tRED */       real_v[24] = real_y[5] + real_y[26];
        /* tRED */       imag_v[24] = imag_y[5] + imag_y[26];
        /* tRED */       real_v[25] -= real_y[26];
        /* tRED */       imag_v[25] -= imag_y[26];
        /* tRED */       real_v[30] = real_y[6];
        /* tRED */       imag_v[30] = imag_y[6];
        /* tRED */       real_v[26] = real_y[6] + real_y[27];
        /* tRED */       imag_v[26] = imag_y[6] + imag_y[27];
        /* tRED */       real_v[30] -= real_y[27];
        /* tRED */       imag_v[30] -= imag_y[27];
        /* tRED */       real_v[27] = real_y[6] + real_y[28];
        /* tRED */       imag_v[27] = imag_y[6] + imag_y[28];
        /* tRED */       real_v[30] -= real_y[28];
        /* tRED */       imag_v[30] -= imag_y[28];
        /* tRED */       real_v[28] = real_y[6] + real_y[29];
        /* tRED */       imag_v[28] = imag_y[6] + imag_y[29];
        /* tRED */       real_v[30] -= real_y[29];
        /* tRED */       imag_v[30] -= imag_y[29];
        /* tRED */       real_v[29] = real_y[6] + real_y[30];
        /* tRED */       imag_v[29] = imag_y[6] + imag_y[30];
        /* tRED */       real_v[30] -= real_y[30];
        /* tRED */       imag_v[30] -= imag_y[30];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       real_y[13] = real_v[13];
        /* tRED */       imag_y[13] = imag_v[13];
        /* tRED */       real_y[14] = real_v[14];
        /* tRED */       imag_y[14] = imag_v[14];
        /* tRED */       real_y[15] = real_v[15];
        /* tRED */       imag_y[15] = imag_v[15];
        /* tRED */       real_y[16] = real_v[16];
        /* tRED */       imag_y[16] = imag_v[16];
        /* tRED */       real_y[17] = real_v[17];
        /* tRED */       imag_y[17] = imag_v[17];
        /* tRED */       real_y[18] = real_v[18];
        /* tRED */       imag_y[18] = imag_v[18];
        /* tRED */       real_y[19] = real_v[19];
        /* tRED */       imag_y[19] = imag_v[19];
        /* tRED */       real_y[20] = real_v[20];
        /* tRED */       imag_y[20] = imag_v[20];
        /* tRED */       real_y[21] = real_v[21];
        /* tRED */       imag_y[21] = imag_v[21];
        /* tRED */       real_y[22] = real_v[22];
        /* tRED */       imag_y[22] = imag_v[22];
        /* tRED */       real_y[23] = real_v[23];
        /* tRED */       imag_y[23] = imag_v[23];
        /* tRED */       real_y[24] = real_v[24];
        /* tRED */       imag_y[24] = imag_v[24];
        /* tRED */       real_y[25] = real_v[25];
        /* tRED */       imag_y[25] = imag_v[25];
        /* tRED */       real_y[26] = real_v[26];
        /* tRED */       imag_y[26] = imag_v[26];
        /* tRED */       real_y[27] = real_v[27];
        /* tRED */       imag_y[27] = imag_v[27];
        /* tRED */       real_y[28] = real_v[28];
        /* tRED */       imag_y[28] = imag_v[28];
        /* tRED */       real_y[29] = real_v[29];
        /* tRED */       imag_y[29] = imag_v[29];
        /* tRED */       real_y[30] = real_v[30];
        /* tRED */       imag_y[30] = imag_v[30];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 3 a 2 c 5  */
        /* tRED */       real_v[11] = real_y[1];
        /* tRED */       imag_v[11] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[11];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[11];
        /* tRED */       real_v[11] -= real_y[11];
        /* tRED */       imag_v[11] -= imag_y[11];
        /* tRED */       real_v[6] = real_y[1] + real_y[16];
        /* tRED */       imag_v[6] = imag_y[1] + imag_y[16];
        /* tRED */       real_v[11] -= real_y[16];
        /* tRED */       imag_v[11] -= imag_y[16];
        /* tRED */       real_v[12] = real_y[2];
        /* tRED */       imag_v[12] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[12];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[12];
        /* tRED */       real_v[12] -= real_y[12];
        /* tRED */       imag_v[12] -= imag_y[12];
        /* tRED */       real_v[7] = real_y[2] + real_y[17];
        /* tRED */       imag_v[7] = imag_y[2] + imag_y[17];
        /* tRED */       real_v[12] -= real_y[17];
        /* tRED */       imag_v[12] -= imag_y[17];
        /* tRED */       real_v[13] = real_y[3];
        /* tRED */       imag_v[13] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[13];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[13];
        /* tRED */       real_v[13] -= real_y[13];
        /* tRED */       imag_v[13] -= imag_y[13];
        /* tRED */       real_v[8] = real_y[3] + real_y[18];
        /* tRED */       imag_v[8] = imag_y[3] + imag_y[18];
        /* tRED */       real_v[13] -= real_y[18];
        /* tRED */       imag_v[13] -= imag_y[18];
        /* tRED */       real_v[14] = real_y[4];
        /* tRED */       imag_v[14] = imag_y[4];
        /* tRED */       real_v[4] = real_y[4] + real_y[14];
        /* tRED */       imag_v[4] = imag_y[4] + imag_y[14];
        /* tRED */       real_v[14] -= real_y[14];
        /* tRED */       imag_v[14] -= imag_y[14];
        /* tRED */       real_v[9] = real_y[4] + real_y[19];
        /* tRED */       imag_v[9] = imag_y[4] + imag_y[19];
        /* tRED */       real_v[14] -= real_y[19];
        /* tRED */       imag_v[14] -= imag_y[19];
        /* tRED */       real_v[15] = real_y[5];
        /* tRED */       imag_v[15] = imag_y[5];
        /* tRED */       real_v[5] = real_y[5] + real_y[15];
        /* tRED */       imag_v[5] = imag_y[5] + imag_y[15];
        /* tRED */       real_v[15] -= real_y[15];
        /* tRED */       imag_v[15] -= imag_y[15];
        /* tRED */       real_v[10] = real_y[5] + real_y[20];
        /* tRED */       imag_v[10] = imag_y[5] + imag_y[20];
        /* tRED */       real_v[15] -= real_y[20];
        /* tRED */       imag_v[15] -= imag_y[20];
        /* tRED */       real_v[26] = real_y[6];
        /* tRED */       imag_v[26] = imag_y[6];
        /* tRED */       real_v[16] = real_y[6] + real_y[21];
        /* tRED */       imag_v[16] = imag_y[6] + imag_y[21];
        /* tRED */       real_v[26] -= real_y[21];
        /* tRED */       imag_v[26] -= imag_y[21];
        /* tRED */       real_v[21] = real_y[6] + real_y[26];
        /* tRED */       imag_v[21] = imag_y[6] + imag_y[26];
        /* tRED */       real_v[26] -= real_y[26];
        /* tRED */       imag_v[26] -= imag_y[26];
        /* tRED */       real_v[27] = real_y[7];
        /* tRED */       imag_v[27] = imag_y[7];
        /* tRED */       real_v[17] = real_y[7] + real_y[22];
        /* tRED */       imag_v[17] = imag_y[7] + imag_y[22];
        /* tRED */       real_v[27] -= real_y[22];
        /* tRED */       imag_v[27] -= imag_y[22];
        /* tRED */       real_v[22] = real_y[7] + real_y[27];
        /* tRED */       imag_v[22] = imag_y[7] + imag_y[27];
        /* tRED */       real_v[27] -= real_y[27];
        /* tRED */       imag_v[27] -= imag_y[27];
        /* tRED */       real_v[28] = real_y[8];
        /* tRED */       imag_v[28] = imag_y[8];
        /* tRED */       real_v[18] = real_y[8] + real_y[23];
        /* tRED */       imag_v[18] = imag_y[8] + imag_y[23];
        /* tRED */       real_v[28] -= real_y[23];
        /* tRED */       imag_v[28] -= imag_y[23];
        /* tRED */       real_v[23] = real_y[8] + real_y[28];
        /* tRED */       imag_v[23] = imag_y[8] + imag_y[28];
        /* tRED */       real_v[28] -= real_y[28];
        /* tRED */       imag_v[28] -= imag_y[28];
        /* tRED */       real_v[29] = real_y[9];
        /* tRED */       imag_v[29] = imag_y[9];
        /* tRED */       real_v[19] = real_y[9] + real_y[24];
        /* tRED */       imag_v[19] = imag_y[9] + imag_y[24];
        /* tRED */       real_v[29] -= real_y[24];
        /* tRED */       imag_v[29] -= imag_y[24];
        /* tRED */       real_v[24] = real_y[9] + real_y[29];
        /* tRED */       imag_v[24] = imag_y[9] + imag_y[29];
        /* tRED */       real_v[29] -= real_y[29];
        /* tRED */       imag_v[29] -= imag_y[29];
        /* tRED */       real_v[30] = real_y[10];
        /* tRED */       imag_v[30] = imag_y[10];
        /* tRED */       real_v[20] = real_y[10] + real_y[25];
        /* tRED */       imag_v[20] = imag_y[10] + imag_y[25];
        /* tRED */       real_v[30] -= real_y[25];
        /* tRED */       imag_v[30] -= imag_y[25];
        /* tRED */       real_v[25] = real_y[10] + real_y[30];
        /* tRED */       imag_v[25] = imag_y[10] + imag_y[30];
        /* tRED */       real_v[30] -= real_y[30];
        /* tRED */       imag_v[30] -= imag_y[30];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       real_y[13] = real_v[13];
        /* tRED */       imag_y[13] = imag_v[13];
        /* tRED */       real_y[14] = real_v[14];
        /* tRED */       imag_y[14] = imag_v[14];
        /* tRED */       real_y[15] = real_v[15];
        /* tRED */       imag_y[15] = imag_v[15];
        /* tRED */       real_y[16] = real_v[16];
        /* tRED */       imag_y[16] = imag_v[16];
        /* tRED */       real_y[17] = real_v[17];
        /* tRED */       imag_y[17] = imag_v[17];
        /* tRED */       real_y[18] = real_v[18];
        /* tRED */       imag_y[18] = imag_v[18];
        /* tRED */       real_y[19] = real_v[19];
        /* tRED */       imag_y[19] = imag_v[19];
        /* tRED */       real_y[20] = real_v[20];
        /* tRED */       imag_y[20] = imag_v[20];
        /* tRED */       real_y[21] = real_v[21];
        /* tRED */       imag_y[21] = imag_v[21];
        /* tRED */       real_y[22] = real_v[22];
        /* tRED */       imag_y[22] = imag_v[22];
        /* tRED */       real_y[23] = real_v[23];
        /* tRED */       imag_y[23] = imag_v[23];
        /* tRED */       real_y[24] = real_v[24];
        /* tRED */       imag_y[24] = imag_v[24];
        /* tRED */       real_y[25] = real_v[25];
        /* tRED */       imag_y[25] = imag_v[25];
        /* tRED */       real_y[26] = real_v[26];
        /* tRED */       imag_y[26] = imag_v[26];
        /* tRED */       real_y[27] = real_v[27];
        /* tRED */       imag_y[27] = imag_v[27];
        /* tRED */       real_y[28] = real_v[28];
        /* tRED */       imag_y[28] = imag_v[28];
        /* tRED */       real_y[29] = real_v[29];
        /* tRED */       imag_y[29] = imag_v[29];
        /* tRED */       real_y[30] = real_v[30];
        /* tRED */       imag_y[30] = imag_v[30];
        /* tRED */       /*  tRED exit  */
        /* tRED */       /*  tRED entry p 2 a 1 c 15  */
        /* tRED */       real_v[16] = real_y[1];
        /* tRED */       imag_v[16] = imag_y[1];
        /* tRED */       real_v[1] = real_y[1] + real_y[16];
        /* tRED */       imag_v[1] = imag_y[1] + imag_y[16];
        /* tRED */       real_v[16] -= real_y[16];
        /* tRED */       imag_v[16] -= imag_y[16];
        /* tRED */       real_v[17] = real_y[2];
        /* tRED */       imag_v[17] = imag_y[2];
        /* tRED */       real_v[2] = real_y[2] + real_y[17];
        /* tRED */       imag_v[2] = imag_y[2] + imag_y[17];
        /* tRED */       real_v[17] -= real_y[17];
        /* tRED */       imag_v[17] -= imag_y[17];
        /* tRED */       real_v[18] = real_y[3];
        /* tRED */       imag_v[18] = imag_y[3];
        /* tRED */       real_v[3] = real_y[3] + real_y[18];
        /* tRED */       imag_v[3] = imag_y[3] + imag_y[18];
        /* tRED */       real_v[18] -= real_y[18];
        /* tRED */       imag_v[18] -= imag_y[18];
        /* tRED */       real_v[19] = real_y[4];
        /* tRED */       imag_v[19] = imag_y[4];
        /* tRED */       real_v[4] = real_y[4] + real_y[19];
        /* tRED */       imag_v[4] = imag_y[4] + imag_y[19];
        /* tRED */       real_v[19] -= real_y[19];
        /* tRED */       imag_v[19] -= imag_y[19];
        /* tRED */       real_v[20] = real_y[5];
        /* tRED */       imag_v[20] = imag_y[5];
        /* tRED */       real_v[5] = real_y[5] + real_y[20];
        /* tRED */       imag_v[5] = imag_y[5] + imag_y[20];
        /* tRED */       real_v[20] -= real_y[20];
        /* tRED */       imag_v[20] -= imag_y[20];
        /* tRED */       real_v[21] = real_y[6];
        /* tRED */       imag_v[21] = imag_y[6];
        /* tRED */       real_v[6] = real_y[6] + real_y[21];
        /* tRED */       imag_v[6] = imag_y[6] + imag_y[21];
        /* tRED */       real_v[21] -= real_y[21];
        /* tRED */       imag_v[21] -= imag_y[21];
        /* tRED */       real_v[22] = real_y[7];
        /* tRED */       imag_v[22] = imag_y[7];
        /* tRED */       real_v[7] = real_y[7] + real_y[22];
        /* tRED */       imag_v[7] = imag_y[7] + imag_y[22];
        /* tRED */       real_v[22] -= real_y[22];
        /* tRED */       imag_v[22] -= imag_y[22];
        /* tRED */       real_v[23] = real_y[8];
        /* tRED */       imag_v[23] = imag_y[8];
        /* tRED */       real_v[8] = real_y[8] + real_y[23];
        /* tRED */       imag_v[8] = imag_y[8] + imag_y[23];
        /* tRED */       real_v[23] -= real_y[23];
        /* tRED */       imag_v[23] -= imag_y[23];
        /* tRED */       real_v[24] = real_y[9];
        /* tRED */       imag_v[24] = imag_y[9];
        /* tRED */       real_v[9] = real_y[9] + real_y[24];
        /* tRED */       imag_v[9] = imag_y[9] + imag_y[24];
        /* tRED */       real_v[24] -= real_y[24];
        /* tRED */       imag_v[24] -= imag_y[24];
        /* tRED */       real_v[25] = real_y[10];
        /* tRED */       imag_v[25] = imag_y[10];
        /* tRED */       real_v[10] = real_y[10] + real_y[25];
        /* tRED */       imag_v[10] = imag_y[10] + imag_y[25];
        /* tRED */       real_v[25] -= real_y[25];
        /* tRED */       imag_v[25] -= imag_y[25];
        /* tRED */       real_v[26] = real_y[11];
        /* tRED */       imag_v[26] = imag_y[11];
        /* tRED */       real_v[11] = real_y[11] + real_y[26];
        /* tRED */       imag_v[11] = imag_y[11] + imag_y[26];
        /* tRED */       real_v[26] -= real_y[26];
        /* tRED */       imag_v[26] -= imag_y[26];
        /* tRED */       real_v[27] = real_y[12];
        /* tRED */       imag_v[27] = imag_y[12];
        /* tRED */       real_v[12] = real_y[12] + real_y[27];
        /* tRED */       imag_v[12] = imag_y[12] + imag_y[27];
        /* tRED */       real_v[27] -= real_y[27];
        /* tRED */       imag_v[27] -= imag_y[27];
        /* tRED */       real_v[28] = real_y[13];
        /* tRED */       imag_v[28] = imag_y[13];
        /* tRED */       real_v[13] = real_y[13] + real_y[28];
        /* tRED */       imag_v[13] = imag_y[13] + imag_y[28];
        /* tRED */       real_v[28] -= real_y[28];
        /* tRED */       imag_v[28] -= imag_y[28];
        /* tRED */       real_v[29] = real_y[14];
        /* tRED */       imag_v[29] = imag_y[14];
        /* tRED */       real_v[14] = real_y[14] + real_y[29];
        /* tRED */       imag_v[14] = imag_y[14] + imag_y[29];
        /* tRED */       real_v[29] -= real_y[29];
        /* tRED */       imag_v[29] -= imag_y[29];
        /* tRED */       real_v[30] = real_y[15];
        /* tRED */       imag_v[30] = imag_y[15];
        /* tRED */       real_v[15] = real_y[15] + real_y[30];
        /* tRED */       imag_v[15] = imag_y[15] + imag_y[30];
        /* tRED */       real_v[30] -= real_y[30];
        /* tRED */       imag_v[30] -= imag_y[30];
        /* tRED */       real_y[1] = real_v[1];
        /* tRED */       imag_y[1] = imag_v[1];
        /* tRED */       real_y[2] = real_v[2];
        /* tRED */       imag_y[2] = imag_v[2];
        /* tRED */       real_y[3] = real_v[3];
        /* tRED */       imag_y[3] = imag_v[3];
        /* tRED */       real_y[4] = real_v[4];
        /* tRED */       imag_y[4] = imag_v[4];
        /* tRED */       real_y[5] = real_v[5];
        /* tRED */       imag_y[5] = imag_v[5];
        /* tRED */       real_y[6] = real_v[6];
        /* tRED */       imag_y[6] = imag_v[6];
        /* tRED */       real_y[7] = real_v[7];
        /* tRED */       imag_y[7] = imag_v[7];
        /* tRED */       real_y[8] = real_v[8];
        /* tRED */       imag_y[8] = imag_v[8];
        /* tRED */       real_y[9] = real_v[9];
        /* tRED */       imag_y[9] = imag_v[9];
        /* tRED */       real_y[10] = real_v[10];
        /* tRED */       imag_y[10] = imag_v[10];
        /* tRED */       real_y[11] = real_v[11];
        /* tRED */       imag_y[11] = imag_v[11];
        /* tRED */       real_y[12] = real_v[12];
        /* tRED */       imag_y[12] = imag_v[12];
        /* tRED */       real_y[13] = real_v[13];
        /* tRED */       imag_y[13] = imag_v[13];
        /* tRED */       real_y[14] = real_v[14];
        /* tRED */       imag_y[14] = imag_v[14];
        /* tRED */       real_y[15] = real_v[15];
        /* tRED */       imag_y[15] = imag_v[15];
        /* tRED */       real_y[16] = real_v[16];
        /* tRED */       imag_y[16] = imag_v[16];
        /* tRED */       real_y[17] = real_v[17];
        /* tRED */       imag_y[17] = imag_v[17];
        /* tRED */       real_y[18] = real_v[18];
        /* tRED */       imag_y[18] = imag_v[18];
        /* tRED */       real_y[19] = real_v[19];
        /* tRED */       imag_y[19] = imag_v[19];
        /* tRED */       real_y[20] = real_v[20];
        /* tRED */       imag_y[20] = imag_v[20];
        /* tRED */       real_y[21] = real_v[21];
        /* tRED */       imag_y[21] = imag_v[21];
        /* tRED */       real_y[22] = real_v[22];
        /* tRED */       imag_y[22] = imag_v[22];
        /* tRED */       real_y[23] = real_v[23];
        /* tRED */       imag_y[23] = imag_v[23];
        /* tRED */       real_y[24] = real_v[24];
        /* tRED */       imag_y[24] = imag_v[24];
        /* tRED */       real_y[25] = real_v[25];
        /* tRED */       imag_y[25] = imag_v[25];
        /* tRED */       real_y[26] = real_v[26];
        /* tRED */       imag_y[26] = imag_v[26];
        /* tRED */       real_y[27] = real_v[27];
        /* tRED */       imag_y[27] = imag_v[27];
        /* tRED */       real_y[28] = real_v[28];
        /* tRED */       imag_y[28] = imag_v[28];
        /* tRED */       real_y[29] = real_v[29];
        /* tRED */       imag_y[29] = imag_v[29];
        /* tRED */       real_y[30] = real_v[30];
        /* tRED */       imag_y[30] = imag_v[30];
        /* tRED */       /*  tRED exit  */
      /* tKRED */       /*  tKRED exit */


        //data[op]		
        for (int px = 0; px < FFTLENGTH; px++)
        {
            real_x[px] = real_y[active_op[px]];
            imag_x[px] = imag_y[active_op[px]];
        }

        for (int px = 0; px < FFTLENGTH; px++)
        {
            real[ind[px]] = real_x[px];
            imag[ind[px]] = imag_x[px];
        }

        //
        //  CRT mapping.
        //
        IncIndices(ind);
    }
}
#undef FFTLENGTH
