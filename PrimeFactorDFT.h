#pragma once
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
#include <vector>

#if 0   // windows
#define WIN
typedef unsigned __int64  u64;
typedef __int64  s64;
#else   //linux
#define NOTWIN
#include <stdlib.h>
#include <cstdint>
typedef  uint64_t  u64;
typedef  int64_t  s64;
#endif

typedef unsigned int uint;
typedef std::vector<uint> factorSeq;
typedef double Data;


class BasicDFT {

public:
	BasicDFT() { count = 0; };
	virtual void Evaluate(Data *read, Data *imag) = 0;

protected:
	std::vector<s64> indices;
	s64 count;

	void IncIndices(std::vector<s64>& ind)
	{
		s64 tmp = ind[ind.size() - 1];
		for (std::size_t i = ind.size() - 1; i > 0; i--)
		{
			ind[i] = ind[i - 1] + 1;
		}
		ind[0] = tmp + 1;
	}


};

class DFT2 : protected BasicDFT {
public:
	DFT2(int  Rotation, s64 Count, std::vector<s64> startIndices)
	{
		count = Count;

		indices = startIndices;

	};

	void Evaluate(Data* real, Data* imag);
private:

};


#undef FFTLENGTH
#define FFTLENGTH 3

class DFT3 : protected BasicDFT {
public:
	DFT3(int  Rotation, s64 Count, std::vector<s64> startIndices) :
		u{
		/*real*/
		-1.500000000000000,
		/* imag */
		0.866025403784439},
		op{ 0, 2, 1},
		ip{ 0, 1, 2} // not used
	{
		int Rotations[FFTLENGTH] = { 0, 1, 2};

		count = Count;
		indices = startIndices;

		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];
	};

	void Evaluate(Data *real, Data *imag);
private:
	const Data  u[2];
	const unsigned int  ip[FFTLENGTH];
	const unsigned int	op[FFTLENGTH];
	unsigned int active_op[FFTLENGTH];

};

#undef FFTLENGTH
#define FFTLENGTH 5

class DFT5 : protected BasicDFT {
public:
	DFT5(int  Rotation, s64 Count, std::vector<s64> startIndices):
		u{ 
		/* real */
		-1.250000000000000,
		-0.559016994374947,
		/* imaginary*/
		-1.538841768587627,
		-0.363271264002681,
		0.951056516295154 
	},
	ip{ 0, 1, 2, 4, 3 },
	op{ 0, 4, 1, 3, 2 }
	{
		int Rotations[FFTLENGTH] = { 0, 1, 2, 3, 4 };

		count = Count;
		indices = startIndices;

		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];
	};

	void Evaluate(Data* real, Data* imag);
private:
	const Data  u[5];
	const unsigned int  ip[FFTLENGTH];
	const unsigned int	op[FFTLENGTH];
	unsigned int active_op[FFTLENGTH];

};

#undef FFTLENGTH
#define FFTLENGTH 7

class DFT7 : protected BasicDFT {
public:
	DFT7(int  Rotation, s64 Count, std::vector<s64> startIndices) :
		u{ 
		/* real */
		-1.166666666666667,
		/* Imag */
		0.440958551844098,
		/* real */
	   -0.678447933946105,
		0.846010735815048,
	   -0.055854267289648,
		/* Imag */
	   -1.408811651299382,
	   -0.193096429713794,
		0.533969360337725 },
		ip{ 0, 1, 4, 2, 6, 3, 5 },
		op{ 0, 6, 5, 1, 4, 2, 3 }
	{
		int Rotations[FFTLENGTH] = { 0, 1, 2, 3, 4, 5, 6 };

		count = Count;
		indices = startIndices;

		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];
	};

	void Evaluate(Data* real, Data* imag);
private:
	const Data  u[8];
	const unsigned int  ip[FFTLENGTH];
	const unsigned int	op[FFTLENGTH];
	unsigned int active_op[FFTLENGTH];

};
#undef FFTLENGTH
#define FFTLENGTH 11
class DFT11 : protected BasicDFT {
public:
	DFT11(int  Rotation, s64 Count, std::vector<s64> startIndices) :
		op{ 0, 10, 1, 8, 7, 9, 4, 2, 3, 6, 5 },
		ip{ 0,  1, 9, 4, 3,	5, 10,2, 7,	8, 6 },
		u{
		 -1.100000000000000,
		0.331662479035540,
		/* pure real*/
		0.253097611605959,
		-1.288200610773679,
		0.304632239669212, 
		-0.391339615511917,
		-2.871022253392850,
		1.374907986616384,
		0.817178135341212,
		1.800746506445679,
		 -0.859492973614498,
		/* pure imaginary */
		 -2.373470454748280,
		 -0.024836393087493,
		 0.474017017512829,
		 0.742183927770612 ,
		 1.406473309094609,
		 -1.191364552195948,
		 0.708088885039503 ,
		 0.258908260614168 ,
		 -0.049929922194110
		}
	{

		int Rotations[FFTLENGTH] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

		count = Count;
		indices = startIndices;

		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];

	}
	

	void Evaluate(Data* real, Data* imag);

private:

	const Data  u[20];
	const unsigned int  ip[FFTLENGTH];
	const unsigned int	op[FFTLENGTH];
	unsigned int active_op[FFTLENGTH];

};

#undef FFTLENGTH
#define FFTLENGTH 13

class DFT13 : protected BasicDFT {

public:
	DFT13(int  Rotation, s64 Count, std::vector<s64> startIndices) :
		op{ 0,12,1,10,5,3,2,8,9,11,4,7,6 },
		ip{ 0,1,3,9,5,2,6,12,10,4,8,11,7 },
		u{
		-1.083333333333333,
		-0.300462606288666,
		-0.749279330626139,
		/* imag */
		0.401002128321867,
		0.174138601152136,
		/* real */
		1.007074065727533,
		0.731245990975348,
		-0.579440018900960,
		0.531932498429674,
		-0.508814921720398,
		-0.007705858903092,
		/* imag  */
		 -2.511393318389568,
		-1.823546408682421,
		1.444979909023996,
		-1.344056915177370,
		-0.975932420775946,
		0.773329778651105,
		1.927725116783469,
		1.399739414729183,
		-1.109154843837551
		} 
	{
		count = Count;
		indices = startIndices;

		int Rotations[FFTLENGTH] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];
	}
	void Evaluate(Data* real, Data* imag);
private:

	const Data  u[20];

	const unsigned int  ip[FFTLENGTH];
	const unsigned int	op[FFTLENGTH];

	unsigned int active_op[FFTLENGTH];

};
#undef FFTLENGTH
#define FFTLENGTH 17

class DFT17 : protected BasicDFT {

public:
	DFT17(int  Rotation, s64 Count, std::vector<s64> startIndices) :
		u{
		/* real */
			-1.062500000000000,
			-0.257694101601104,
			0.723407977286057,
			-0.089055591620606,
			-0.317176192832725,
			0.924380996081242,
			0.676798496730885,
			-0.440889073481754,
			-1.517002366671939,
			-0.797601020823318,
			1.281092943422807,
			0.296310685295348,
			0.060401262046216,
			-0.420101934970527,
			/* imag */
			1.462686052158509,
			 2.709842506062867,
			 -1.124438635937869,
			 -1.808356521480244,
			 2.958485673330231 ,
			 0.222952651355246 ,
			 -0.906077574510765 ,
			 -2.491481444635763 ,
			 0.634492510107882 ,
			 2.681907643666417 ,
			 -1.890642422994411 ,
			 0.499530681019060 ,
			 0.524082025323152 ,
			 -1.205277132872840 ,
			 0.867029716652200 ,
			 0.032526452324592 ,
			 1.423638194299944 ,
			 -1.356975842482187 ,
			 -2.072296847912463 ,
			 -0.409600041534228 ,
			 0.312453977459405 ,
			 0.642137248078546 ,
			 -0.876604270228695 ,
			 -0.544991184003723 ,
			 0.436775561093086 ,
			 0.533921625167910 ,
			 0.361241666187153 

	},
		ip{ 0,1,3,9,10,13,5,15,11,16,14,8,7,4,12,2,6 },
		op{ 0,16,14,1,12,5,15,11,10,2,3,7,13,4,9,6,8 }
	{
		count = Count;
		indices = startIndices;

		int Rotations[FFTLENGTH] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];
	}
	void Evaluate(Data* real, Data* imag);

private:

	const Data u[41];
	const unsigned int  ip[FFTLENGTH];
	const unsigned int	op[FFTLENGTH];
	unsigned int active_op[FFTLENGTH];

};


#undef FFTLENGTH
#define FFTLENGTH 19

class DFT19 : protected BasicDFT {

public:
	DFT19(int  Rotation, s64 Count, std::vector<s64> startIndices) :
		u{
			/* real */
			-1.055555555555556,
			/* imag */
			 0.242161052418926,
			/* real */
			0.798693520987127,
			0.177211053261099,
			-0.325301524749409,
			/* imag */
			-0.834854293606883,
			-0.488430732011460 ,
			0.441095008539447,
			/* real */
			0.435557826755211, 
			0.231321070206014, 
			-0.421744310987423,
			-0.002942234699835,
			0.822164874728519, 
			-1.524433450111895,
			-0.208976399520928,
			0.861151390434984, 
			0.242730955603080, 
			-3.304023490019812,
			0.362958541118895, 
			-0.007448223561695,
			-0.146469026482520,
			-0.079929573634415,
			0.827286205097097, 
			/* imag */
			 0.490936114006330,
			 0.364666773063770,
			 -0.318086136404994,
			 -0.314562985092245,
			 5.737605861191472,
			 1.814966178076640,
			 -0.151292354260379,
			 0.313346892339940,
			 -0.143269744887002,
			 2.935067055771071,
			 -0.768634097360990,
			 -0.071124806267797,
			 0.001579748021685,
			 0.152610909993082,
			 -2.890890972320848
	},
	ip{ 0,1,17,4,11,16,	6,7,5,9,18,2,15,8,	3,13,12,14,10 },
	op{ 0,18,1,4,11,16,	14,	15,	3,17,8,	12,	6,5,7,2,13,	10,	9 }
	{
		count = Count;
		indices = startIndices;

		int Rotations[FFTLENGTH] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];
	}

	void Evaluate(Data* real, Data* imag);
private:
	const Data u[39];
	const unsigned int  ip[FFTLENGTH];
	const unsigned int	op[FFTLENGTH];
	unsigned int active_op[FFTLENGTH];

};

#undef FFTLENGTH
#define FFTLENGTH 31


class DFT31 : protected BasicDFT {

public:
	DFT31(int  Rotation, s64 Count, std::vector<s64> startIndices) :
		op{ 0,30,29,1,28,25,5,18,27,22,24,8,4,6,17,11,26,2,21,19,23,9,7,12,3,20,10,13,16,14,15 },
		ip{ 0,1,16,8,4,2,25,28,14,7,19,5,18,9,20,10,30,15,23,27,29,6,3,17,24,12,26,13,22,11,21 },
		u{ 
		/* real */
		-1.033333333333333,	/*  0 */
		/* imag */
		0.185592145427667 ,	/*  1 */
		/* real */
		 0.251026872929094, 	/*  2 */
		 0.638094290379888, 	/*  3 */
		 -0.296373721102994,	/*  4 */
		/* imag */
		 -0.462201919825109,	/*  5 */
		 0.155909426230360 ,	/*  6 */
		 0.102097497864916 ,	/*  7 */
		/* real */
		 -0.100498239164838, 	/*  8 */
		 -0.217421331841463, 	/*  9 */
		 -0.325082164955763, 	/*  10 */
		 0.798589508696894, 	/*  11 */
		 -0.780994042074251,	/*  12 */
		  -0.256086011899669,	/*  13 */
		 0.169494392220932,	/*  14 */
		 0.711997889018157, 	/*  15 */
		 -0.060064820876732,	/*  16 */
		/* imag */
		 -1.235197570427205 ,	/*  17 */
		 -0.271691369288525 ,	/*  18 */
		 0.541789612349592 ,	/*  19 */
		 0.329410560797314 ,	/*  20 */
		 1.317497505049809 ,	/*  21 */
		 -0.599508803858381 ,	/*  22 */
		 0.093899154219231 ,	/*  23 */
		 -0.176199088841836 ,	/*  24 */
		 0.028003825226279 ,	/*  25 */
		/* real */
		 1.316699050305790, 	/*  26 */
		 1.330315270540553, 	/*  27 */
		 -0.385122753006171,	/*  28 */
		 -2.958666546021397,	/*  29 */
		 -2.535301995146201,	/*  30 */
		 2.013474028487015, 	/*  31 */
		 1.081897731187396, 	/*  32 */
		 0.136705213653014, 	/*  33 */
		 -0.569390844064251,	/*  34 */
		 -0.262247009112805,	/*  35 */
		 2.009855570455675, 	/*  36 */
		 -1.159348599757857,	/*  37 */
		 0.629367699727360, 	/*  38 */
		 1.229312102919654, 	/*  39 */
		 -1.479874670425178,	/*  40 */
		 -0.058279061554516,	/*  41 */
		 -0.908786032252333,	/*  42 */
		 0.721257672797977, 	/*  43 */
		 -0.351484013730995,	/*  44 */
		 -1.113390280332076,	/*  45 */
		 0.514823784254676, 	/*  46 */
		 0.776432948764679, 	/*  47 */
		 0.435329964075516, 	/*  48 */
		 -0.177866452687279,	/*  49 */
		 -0.341206223210960,	/*  50 */
		 0.257360272866440, 	/*  51 */
		 -0.050622276244575,	/*  52 */
		/* imag */
		 -2.745673340229639 ,	/*  53 */
		 2.685177424507523 ,	/*  54 */
		 0.880463026400118 ,	/*  55 */
		 -5.028851220636894 ,	/*  56 */
		 -0.345528375980267 ,	/*  57 */
		 1.463210769729252 ,	/*  58 */
		 3.328421083558774 ,	/*  59 */
		 -0.237219367348867 ,	/*  60 */
		 -1.086975102467855 ,	/*  61 */
		 -1.665522956385442 ,	/*  62 */
		 1.628826188810638 ,	/*  63 */
		 0.534088072762272 ,	/*  64 */
		 -3.050496586573981 ,	/*  65 */
		 -0.209597199290132 ,	/*  66 */
		 0.887582325001072 ,	/*  67 */
		 2.019017208624242 ,	/*  68 */
		 -0.143897052948668 ,	/*  69 */
		 -0.659358110687783 ,	/*  70 */
		 1.470398765538361 ,	/*  71 */
		 -1.438001204439387 ,	/*  72 */
		 -0.471517033054130 ,	/*  73 */
		 2.693115935736959 ,	/*  74 */
		 0.185041858423467 ,	/*  75 */
		 -0.783597698243441 ,	/*  76 */
		 -1.782479430727672 ,	/*  77 */
		 0.127038806765845 ,	/*  78 */
		 0.582111071051880 	/*  79 */
		}


	{
		count = Count;
		indices = startIndices;
#ifdef WIN
		_Rotation = std::abs(Rotation);
#endif
#ifdef NOTWIN
		_Rotation = abs(Rotation);
#endif

		int Rotations[FFTLENGTH] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
									 17,18,19,20,21,22,23,24,25,26,27,28,29,30 };
		for (int i = 0; i < FFTLENGTH; i++)
		{
			Rotations[i] *= Rotation;
		}
		for (int i = 0; i < FFTLENGTH; i++)
		{
			while (Rotations[i] < 0) Rotations[i] += FFTLENGTH;
			while (Rotations[i] >= FFTLENGTH) Rotations[i] -= FFTLENGTH;
		}
		for (int i = 0; i < FFTLENGTH; i++)
			active_op[i] = op[Rotations[i]];


	}
	void Evaluate(Data *real, Data* imag);

private:

	const Data  u[80];
	const unsigned int  ip[31];
	const unsigned int	op[31];
	unsigned int active_op[31];
	unsigned int _Rotation;

};


class PrimeFactorDFT
{
public:
	
	PrimeFactorDFT() {};

	void SetFactors(factorSeq& _factors) {
		factors = _factors;
		state = ValidateFactors(factors);
		CleanUpDFT(DFTs);
		if (state > 0) {
			InitRotations();
			InitDFT(factors, DFTs);
		}
	};

	void GetFactors(factorSeq& _factors) {_factors = factors;};

	int CalcFactors(uint length, factorSeq& _factors, int factorCount = 0);

	/*
	*  Based of the factors provided.
	*  if > 0 the length of the FFT.
	*  if == 0 no factors provided.
	*  if == -1 invalid/unsupported factors provided.
	*  if == -2 duplicated factor  provided.
	*/
	s64 Status() { return state; };

	void forwardFFT(Data* real, Data *imag);
	void InverseFFT(Data* real, Data *imag);
	void ScaledInverseFFT(Data* real, Data *imag);

private:

	s64 ValidateFactors(factorSeq& _factors);
	s64 state;
	void InitDFT(factorSeq& _factors, std::vector<BasicDFT*> &_DTFs);
	void CleanUpDFT(std::vector<BasicDFT*> &_DTFs);
	void InitRotations();
	void InitIndices(std::vector<s64>& indices, int fftlength, s64 length);
	factorSeq factors;
	std::vector<int>  Rotations;
	std::vector<BasicDFT*> DFTs;
};

