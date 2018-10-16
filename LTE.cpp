#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Rate 3        // coderate^-1
#define termi 3       // number of termination systematic bits
#define N 2051        // interleaver length + termi

double **alpha, **beta, **gamma, variance, nSymbol0[16], nSymbol1[16], bSymbol0[16], bSymbol1[16];
char nextState0[8], nextState1[8], preState0[8], preState1[8];


/* ***************************************************************************** */
/* Copyright:      Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*                 Makoto Matsumoto, Hiroshima University                        */
/* Notice:         This code can be used freely for personal, academic,          */
/*                 or non-commercial purposes. For commercial purposes,          */
/*                 please contact P. L'Ecuyer at: lecuyer@iro.umontreal.ca       */
/* ***************************************************************************** */

#define W 32
#define R 32
#define M1 3
#define M2 24
#define M3 10

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define Identity(v) (v)

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000001fU]
#define VM2           STATE[(state_i+M2) & 0x0000001fU]
#define VM3           STATE[(state_i+M3) & 0x0000001fU]
#define VRm1          STATE[(state_i+31) & 0x0000001fU]
#define newV0         STATE[(state_i+31) & 0x0000001fU]
#define newV1         STATE[state_i                   ]

#define FACT 2.32830643653869628906e-10  // einai o ari8mos 2^-32

static unsigned state_i = 0;
static unsigned STATE[R];
static unsigned z0, z1, z2;

void InitWELLRNG1024a(unsigned *init)
{
	int j;
	state_i = 0;
	for (j = 0; j<R; j++)
		STATE[j] = init[j];
}

int WELLRNG1024a(void)
{
	// the period of the random number generator is 2^1024-1
	z0 = VRm1;
	z1 = Identity(V0) ^ MAT0POS(8, VM1);
	z2 = MAT0NEG(-19, VM2) ^ MAT0NEG(-14, VM3);
	newV1 = z1^z2;
	newV0 = MAT0NEG(-11, z0) ^ MAT0NEG(-7, z1) ^ MAT0NEG(-13, z2);
	state_i = (state_i + 31) & 0x0000001fU;
	return STATE[state_i] - 2147483648;   // epistrefei enan omoiomorfa tyxaio int [-2^31,2^31-1]
}

void zigSet(int *k, double *w, double *f)
{
	double m = 2147483648.0, d = 3.6541528853610088, t = d, v = .00492867323399, q;
	int i;

	q = v / exp(-.5*d*d);
	k[0] = (int)(d / q*m);     k[1] = 0;
	w[0] = q / m;              w[255] = d / m;
	f[0] = 1.;               f[255] = exp(-.5*d*d);

	for (i = 254; i >= 1; i--)
	{
		d = sqrt(-2.*log(v / d + exp(-.5*d*d)));
		k[i + 1] = (int)(d / t*m);
		t = d;
		f[i] = exp(-.5*d*d);
		w[i] = d / m;
	}
}

void mapper(char *help, double *symbol)
{
	char i;
	for (i = 0; i < 8; i++)
		switch (help[i])
		{
			case 0:
				symbol[2 * i] = -1.; symbol[2 * i + 1] = -1.;
			   break;
			case 1:
				symbol[2 * i] = -1.; symbol[2 * i + 1] = 1.;
			   break;
			case 2:
				symbol[2 * i] = 1.; symbol[2 * i + 1] = -1.;
			   break;
			default:
				symbol[2 * i] = 1.; symbol[2 * i + 1] = 1.;
		}

	for (i = 0; i < 16; i++)
		symbol[i] /= variance;
}

void trellis()
{
	char c, i, help[8];

	FILE *fp;
	fp = fopen("trellis.txt", "r");
	if (fp == NULL)
	{
		printf("Can't open trellis.txt.\n");
		system("pause");
		exit(0);
	}

	// next states
	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", nextState0 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", nextState1 + i);

	//fwd outputs
	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, nSymbol0);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, nSymbol1);

	// previous states
	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", preState0 + i);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", preState1 + i);

	//bwd outputs
	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, bSymbol0);

	do { c = getc(fp); } while (c != ':');
	for (i = 0; i < 8; i++)
		fscanf(fp, "%hhd", help + i);
	mapper(help, bSymbol1);

	fclose(fp);
}

void createInt(int *Int, int *deInt)
{
    // random interleaver
	int i, index, *A;

	A = new int[N - termi];
	for (i = 0; i < N - termi; i++)
		A[i] = i;

	for (i = 0; i < N - termi; i++)
	{
	  here:
		index = abs(WELLRNG1024a()) % (N - termi);
		if (A[index] == -1)
			goto here;
		Int[i] = index;
		deInt[index] = i;
		A[index] = -1;
	}
	delete[] A;
}

void source(char *bitStream)
{
	for (int i = 0; i < N - termi; i++)
		bitStream[i] = WELLRNG1024a() > 0 ? 1 : 0;
}

void encode(char *bitStream, double *IntSystematic, char *parity1, char *parity2, int *deInt)
{
	char input, D1[3] = { 0 }, D2[3] = { 0 };
	int i;

	for (i = 0; i < N - termi; i++)
	{
		parity1[i] = bitStream[i] ^ D1[0] ^ D1[1];

		input = bitStream[i] ^ D1[1] ^ D1[2];
		D1[2] = D1[1];
		D1[1] = D1[0];
		D1[0] = input;

		parity2[i] = bitStream[deInt[i]] ^ D2[0] ^ D2[1];

		input = bitStream[deInt[i]] ^ D2[1] ^ D2[2];
		D2[2] = D2[1];
		D2[1] = D2[0];
		D2[0] = input;
	}

	for (; i < N; i++)
	{
		bitStream[i] = D1[1] ^ D1[2];
		parity1[i] = D1[0] ^ D1[2];

		D1[2] = D1[1];
		D1[1] = D1[0];
		D1[0] = 0;

		IntSystematic[i] = D2[1] ^ D2[2] == 0 ? -1. : 1.; // LTE sends the systematic tail bits of the 2nd constituent encoder too	
		parity2[i] = D2[0] ^ D2[2];

		D2[2] = D2[1];
		D2[1] = D2[0];
		D2[0] = 0;
	}
}

void modulate(char *bitStream, char *parity1, char *parity2)
{
	// BPSK modulation: 1 -> +1, 0 -> -1
	for (int i = 0; i < N; i++)
	{
		if (bitStream[i] == 0)
			bitStream[i] = -1;
		if (parity1[i] == 0)
			parity1[i] = -1;
		if (parity2[i] == 0)
			parity2[i] = -1;
	}
}

void wgn(double *GN, int length, int *k, double *w, double *f)
{
	// creates WGN samples by applying the Ziggurat algorithm of Marsaglia & Tsang (2000)
	int randomInt;
	double uni, uni2, x, y, r = 3.6541528853610088;
	short int i;
	int j = 0;

	while (j < length)
	{
		randomInt = WELLRNG1024a();  // gennaei enan tyxaio int, apo ton opoio 8a prokypsei to deigma 8oryvou
		i = WELLRNG1024a() & 255;      // epilegei tyxaia ena apo ta 256 strwmata tou ziggurat, gennwntas enan kainourio int*,
									   // *symfwna me to "An Improved Ziggurat Method to Generate Normal Random Samples" tou Doornik (2005)

		if (abs(randomInt)<k[i])    // to 99.33% twn deigmatwn proerxontai apo auto to if
		{
			GN[j++] = randomInt*w[i];
			continue;
		}

		if (i == 0)
		{
			do
			{
				uni = .5 + WELLRNG1024a()*FACT;
				if (uni == 0) uni = .5 + WELLRNG1024a()*FACT;
				x = -log(uni)*0.27366123732975827;   // o antistrofos tou r
				uni = .5 + WELLRNG1024a()*FACT;
				if (uni == 0) uni = .5 + WELLRNG1024a()*FACT;
				y = -log(uni);
			} while (y + y < x*x);

			GN[j++] = randomInt > 0 ? r + x : -r - x;
			continue;
		}

		uni = randomInt*w[i];
		uni2 = .5 + WELLRNG1024a()*FACT;
		if (f[i] + uni2*(f[i - 1] - f[i]) < exp(-.5*uni*uni))
			GN[j++] = uni;
	}
}

void sisoDec(double *systematic, double *parity, double *La, double *LLR)
{
	// log-APP algorithm
	char numOfStates = 8, s, step = 8;
	int i;
	double x, y;


	//-------------------------- gamma branch metrics calculation ----------------------------//

	for (i = 0; i < N - termi; i++)
	{
		for (s = 0; s < numOfStates; s += step)
		{
			gamma[i][s] = -La[i] / 2 + systematic[i] * nSymbol0[s << 1] + parity[i] * nSymbol0[s << 1 | 1];
			gamma[i][s + 8] = La[i] / 2 + systematic[i] * nSymbol1[s << 1] + parity[i] * nSymbol1[s << 1 | 1];
		}
		if (step != 1) step >>= 1;
	}

	for (; i < N; i++)
	{
		numOfStates >>= 1;

		for (s = 0; s < numOfStates; s++)
		{
			gamma[i][preState0[s]] = systematic[i] * bSymbol0[s << 1] + parity[i] * bSymbol0[s << 1 | 1];
			gamma[i][preState1[s] + 8] = systematic[i] * bSymbol1[s << 1] + parity[i] * bSymbol1[s << 1 | 1];
		}
	}

	//-------------------------- alpha node metrics calculation ----------------------------//

	step = numOfStates = 8;
	for (i = 0; i < termi; i++)
	{
		for (s = 0; s < numOfStates; s += step)
		{
			alpha[i + 1][nextState0[s]] = gamma[i][s] + alpha[i][s];
			alpha[i + 1][nextState1[s]] = gamma[i][s + 8] + alpha[i][s];
		}
		step >>= 1;
	}

	for (; i < N - termi - 1; i++)
	{
		for (s = 0; s < numOfStates; s++)
		{
			x = gamma[i][preState0[s]] + alpha[i][preState0[s]];
			y = gamma[i][preState1[s] + 8] + alpha[i][preState1[s]];

			alpha[i + 1][s] = x > y ? x + log(1 + exp(y - x)) : y + log(1 + exp(x - y)); // max*

			/*
			alpha[i + 1][s] = x > y ? x : y;
			// lookup table
			x = abs(x-y);
			if( x < 0.375 )
			alpha[i+1][s] += 0.6;
			else if( x < 0.75 )
			alpha[i+1][s] += 0.5;
			else if( x < 1. )
			alpha[i+1][s] += 0.4;
			else if( x < 1.5 )
			alpha[i+1][s] += 0.3;
			else if( x < 2.25 )
			alpha[i+1][s] += 0.2;
			else if( x < 3. )
			alpha[i+1][s] += 0.1;
			*/
		}
	}
	numOfStates = 1;

	//-------------------------- beta node metrics calculation -----------------------------//

	for (i = N - 1; i >= N - termi; i--)
	{
		for (s = 0; s<numOfStates; s++)
		{
			beta[i - 1][preState0[s]] = gamma[i][preState0[s]] + beta[i][s];
			beta[i - 1][preState1[s]] = gamma[i][preState1[s] + 8] + beta[i][s];
		}
		numOfStates <<= 1;
	}

	for (; i > 0; i--)
	{
		if (i < termi) step <<= 1;

		for (s = 0; s<numOfStates; s += step)
		{
			x = gamma[i][s] + beta[i][nextState0[s]];
			y = gamma[i][s + 8] + beta[i][nextState1[s]];

			beta[i - 1][s] = x > y ? x + log(1 + exp(y - x)) : y + log(1 + exp(x - y)); // max*

			/*
			beta[i - 1][s] = x > y ? x : y;
			// lookup table
			x = abs(x-y);
			if( x < 0.375 )
			beta[i-1][s] += 0.6;
			else if( x < 0.75 )
			beta[i-1][s] += 0.5;
			else if( x < 1. )
			beta[i-1][s] += 0.4;
			else if( x < 1.5 )
			beta[i-1][s] += 0.3;
			else if( x < 2.25 )
			beta[i-1][s] += 0.2;
			else if( x < 3. )
			beta[i-1][s] += 0.1;
			*/
		}
	}
	step = 8;

	//-------------------------- LLR calculation ----------------------------//

	LLR[0] = beta[0][4] + gamma[0][8] - beta[0][0] - gamma[0][0];

	for (i = 1; i < N - termi; i++)
	{
		if (i <= termi) step >>= 1;

		x = 1.;
		y = alpha[i][0] + gamma[i][8] + beta[i][nextState1[0]];
		for (s = step; s<numOfStates; s += step)
			x += exp(alpha[i][s] + gamma[i][s + 8] + beta[i][nextState1[s]] - y);
		LLR[i] = y + log(x);

		x = 1.;
		y = alpha[i][0] + gamma[i][0] + beta[i][nextState0[0]];
		for (s = step; s<numOfStates; s += step)
			x += exp(alpha[i][s] + gamma[i][s] + beta[i][nextState0[s]] - y);
		LLR[i] -= y + log(x);
	}
}

int main()
{
	bool test = true, genie = true;
	double *GN, *Systematic, *IntSystematic, *noisyParity1, *noisyParity2,
	       *La1, *La2, *LLR,
	       sigma, EsNo, EbNo,
	       wNor[256], fNor[256];
	unsigned seed[R];
	int kNor[256], i, trials = 0, length, bitErrors = 0, frameErrors = 0, maxIters = 8;	
	int *Int, *deInt;
	char *bitStream, *parity1, *parity2, iters, hd;

	//----------------------------------------------------------------------------//

	alpha = new double*[N];
	for (i = 0; i < N; i++)
		alpha[i] = new double[8];
	alpha[0][0] = 0.;

	beta = new double*[N];
	for (i = 0; i < N; i++)
		beta[i] = new double[8];
	beta[N - 1][0] = 0.;

	gamma = new double*[N];
	for (i = 0; i < N; i++)
		gamma[i] = new double[16];

	bitStream = new char[N];
	parity1 = new char[N];
	parity2 = new char[N];

	Int = new int[N - termi];
	deInt = new int[N - termi];

	La1 = new double[N - termi];
	La2 = new double[N - termi];
	LLR = new double[N - termi];

	Systematic = new double[N];
	IntSystematic = new double[N];
	noisyParity1 = new double[N];
	noisyParity2 = new double[N];

	length = Rate*N + termi; // codeword length (+ termi samples for the systematic tail bits of the 2nd constituent encoder)
	GN = new double[length];	

	//----------------------------------------------------------------------------//

	// random number generator initialization
	for (i = 0; i < R; i++)
		seed[i] = 12345 + i;
	InitWELLRNG1024a(seed);
	zigSet(kNor, wNor, fNor);

	EbNo = 1.0; // Eb/No in dB
	EsNo = EbNo + 10 * log10(double((N - termi))/double(length));
	variance = pow(10, -EsNo / 10) / 2.;
	sigma = sqrt(variance);
	trellis();

	//----------------------------------------------------------------------------//

	while (bitErrors < 1000)
	{			
			if ((trials & 31) == 0)
				createInt(Int, deInt); // creates a different random interleaver every 32 trials
			source(bitStream);
			encode(bitStream, IntSystematic, parity1, parity2, deInt);
			modulate(bitStream, parity1, parity2);
			wgn(GN, length, kNor, wNor, fNor);

			for (i = 0; i < N; i++)
			{
				Systematic[i] = bitStream[i] + sigma*GN[i];
				noisyParity1[i] = parity1[i] + sigma*GN[N + i];
				noisyParity2[i] = parity2[i] + sigma*GN[(N << 1) + i];
			}

			for (i = 0; i < N - termi; i++) IntSystematic[Int[i]] = Systematic[i];

			// LTE sends the systematic tail bits of the 2nd constituent encoder too
			for (; i < N; i++) IntSystematic[i] += sigma*GN[(N << 1) + i + termi];	


			// turbo decoder
			for (i = 0; i < N - termi; i++) La1[i] = 0.;
			for (iters = 0; iters < maxIters; iters++)
			{
				sisoDec(Systematic, noisyParity1, La1, LLR);
				for (i = 0; i < N - termi; i++) La2[Int[i]] = LLR[i] - 2 * Systematic[i] / variance - La1[i];

				sisoDec(IntSystematic, noisyParity2, La2, LLR);

				if (genie)
				{
					if (iters > 1 && iters < maxIters - 1)
					{
						for (i = 0; i < N - termi; i++)
						{
							hd = LLR[Int[i]] < 0. ? -1 : 1;
							if (hd != bitStream[i])
								break;
						}
						if (i == N - termi) goto there;
					}
				}

				if (iters < maxIters - 1)
					for (i = 0; i < N - termi; i++) La1[i] = LLR[Int[i]] - 2 * IntSystematic[Int[i]] / variance - La2[Int[i]];
			}

			// bitErrors counting
			for (i = 0; i < N - termi; i++)
			{
				hd = LLR[Int[i]] < 0. ? -1 : 1;
				if (hd != bitStream[i])
				{
					++bitErrors;
					if (test) { frameErrors++; test = false; }
				}
			}
			test = true;

		there:
			trials++;

			printf("%2d %d %d %d %.2e %.2e\n", iters, trials, bitErrors, frameErrors, double(bitErrors) / double(N - termi) / double(trials), double(frameErrors) / double(trials));
	}

	printf("\nEb/No=%.2f  bits=%.4e  BER=%.2e  FER=%.2e\n\n", EbNo, double((N - termi)*trials), double(bitErrors) / double(N - termi) / double(trials), double(frameErrors) / double(trials));

	system("pause");
	return 0;
}
