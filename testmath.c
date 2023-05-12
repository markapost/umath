#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "../umath/fixmath.h"

#define MATRIX_SIZE 3
#define MULTIPLIER 10.0
#define TEST_LENGTH 1000

double dmult(double a, double b)
{
	return a*b;
}

double ddiv(double a, double b)
{
	return a/b;
}

double dsqrt(double a)
{
	return sqrt(a);
}

double dinv(double a)
{
	return 1/a;
}

double dsin(double a)
{
	return sin(a);
}

double dcos(double a)
{
	return cos(a);
}

void PrintMatrix(fix *matrix, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
	{
		printf("[");
		for(j = 0; j < cols; j++)
			printf("%f ", FIX_TO_DOUBLE(matrix[i*rows+j]));
		printf("]\n");
	}
}

int main(int argc, char** argv)
{
	int i, j;
	double a, b, dt, fdt;
	struct timeval timeStart, timeEnd;
	fix fa, fb;
	fix va[MATRIX_SIZE], vb[MATRIX_SIZE], vc[MATRIX_SIZE];
	fix ma[MATRIX_SIZE][MATRIX_SIZE], mb[MATRIX_SIZE*MATRIX_SIZE], mc[MATRIX_SIZE][MATRIX_SIZE],
 		mt[MATRIX_SIZE][MATRIX_SIZE] = {{INT_TO_FIX(1),	INT_TO_FIX(2),	INT_TO_FIX(3)},
						{INT_TO_FIX(4),	INT_TO_FIX(5),	INT_TO_FIX(6)},
						{INT_TO_FIX(7),	INT_TO_FIX(8),	INT_TO_FIX(9)}	};
	double dlist[TEST_LENGTH];
	fix flist[TEST_LENGTH];

	//build a vector of test values
	for(i=0; i<TEST_LENGTH; i++)
	{
		dlist[i] = log2(((double)i+1)*M_PI);
		flist[i] = DOUBLE_TO_FIX(dlist[i]);
	}
	timeStart.tv_sec = time(NULL); //get time in seconds since the Epoch
	srand(timeStart.tv_sec); //initialize random seed

	printf("Byte Sizes: short=%lu int=%lu long=%lu fix=%lu fixrad=%lu float=%lu double=%lu\n",
		(unsigned long)sizeof(short), (unsigned long)sizeof(int), (unsigned long)sizeof(long), (unsigned long)sizeof(fix), (unsigned long)sizeof(fixrad), (unsigned long)sizeof(float), (unsigned long)sizeof(double));

	printf("Fractional bits=%d, scaling=%d (0x%x), step=%f, round by %d (0x%x)\n",
		FIX_FRAC_BITS, FIX_SCALE, FIX_SCALE, FIX_STEP, 1<<(FIX_FRAC_BITS-1), 1<<(FIX_FRAC_BITS-1));

	printf("Max whole %12d (0x%x), Max fraction %12.12lg (0x%x), Maximum=%24.17lg (0x%x)\n",
		FIX_MAX_WHOLE, DOUBLE_TO_FIX(FIX_MAX_WHOLE),
		FIX_MAX_FRACTION, DOUBLE_TO_FIX(FIX_MAX_FRACTION),
		FIX_MAX, DOUBLE_TO_FIX(FIX_MAX));

	printf("Byte Size for Radian Angles=%lu, maximum value %d (0x%x)\n",
		(unsigned long)sizeof(fixrad), FIX_RAD_MAX, FIX_RAD_MAX);

	printf("Constants: FIX_ONE=%g, FIX_PI=%g, FIX_2PI=%g, FIX_E=%g, FIX_ROOT2=%g, FIX_ROOT3=%g, FIX_GOLDEN=%g\n",
		FIX_TO_DOUBLE(FIX_ONE), FIX_TO_DOUBLE(FIX_PI), FIX_TO_DOUBLE(FIX_2PI), FIX_TO_DOUBLE(FIX_E), FIX_TO_DOUBLE(FIX_ROOT2), FIX_TO_DOUBLE(FIX_ROOT3), FIX_TO_DOUBLE(FIX_GOLDEN));

	printf("\nTesting Addition and Subtraction\n");
	for(i = 0; i < 11; i++)
	{
		b = ((double)rand()/(double)RAND_MAX)*exp(i);
		a = ((double)rand()/(double)RAND_MAX)*exp(i);
		fa = DOUBLE_TO_FIX(a);
		fb = DOUBLE_TO_FIX(b);

		//printf("%12f\t+%12f\t=%12f(double)\t%12f(fix)\terr: %.10f%%\n",
		printf("%12f,%12f,%12f,%12f,%.10f,",
			a, b, a+b, FIX_TO_DOUBLE(fa+fb), 100.0*fixabs((a+b)-(FIX_TO_DOUBLE(fa+fb)))/fixabs(a+b));

		//printf("%12f\t-%12f\t=%12f(double)\t%12f(fix)\terr: %.10f%%\n",
		printf("%12f,%12f,%12f,%12f,%.10f\n",
			a, b, a-b, FIX_TO_DOUBLE(fa-fb), 100.0*fixabs((a-b)-(FIX_TO_DOUBLE(fa-fb)))/fixabs(a-b));

	}

	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += dlist[j] + dlist[i];
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nAddition with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += flist[j] + flist[i];
		}
	gettimeofday(&timeEnd, NULL);
	fdt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nAddition with fix: %12f s (%f), improvement %f\n", fdt, FIX_TO_DOUBLE(fa), dt/fdt);
	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += dlist[j] - dlist[i];
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nSubtraction with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += flist[j] - flist[i];
		}
	gettimeofday(&timeEnd, NULL);
	fdt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nSubtraction with fix: %12f s (%f), improvement %f\n", fdt, FIX_TO_DOUBLE(fa), dt/fdt);

	printf("\nTesting Multiplication and Division\n");
	for(i = 0; i < 6; i++)
	{
		b = ((double)rand()/(double)RAND_MAX)*exp(i);
		a = ((double)rand()/(double)RAND_MAX)*exp(i);
		fa = DOUBLE_TO_FIX(a);
		fb = DOUBLE_TO_FIX(b);

		//printf("%12f\t*%12f\t=%12f(double)\t%12f(fix)\terr: %.10f%%\n",
		printf("%12f,%12f,%12f,%12f,%.10f,",
			a, b, a*b, FIX_TO_DOUBLE(fixmult(fa,fb)), 100.0*fixabs((a*b)-(FIX_TO_DOUBLE(fixmult(fa,fb))))/fixabs(a*b));

		//printf("%12f\t/%12f\t=%12f(double)\t%12f(fix)\terr: %.10f%%\n",
		printf("%12f,%12f,%12f,%12f,%.10f\n",
			a, b, a/b, FIX_TO_DOUBLE(fixdiv(fa,fb)), 100.0*fixabs((a/b)-(FIX_TO_DOUBLE(fixdiv(fa,fb))))/fixabs(a/b));
	}

	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += dmult(dlist[j], dlist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nMultiplication with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += fixmult(flist[j],flist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	fdt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nMultiplication with fix: %12f s (%f), improvement %f\n", fdt, FIX_TO_DOUBLE(fa), dt/fdt);
	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += ddiv(dlist[j], dlist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nDivision with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += fixdiv(flist[j],flist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	fdt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nDivision with fix: %12f s (%f), improvement %f\n", fdt, FIX_TO_DOUBLE(fa), dt/fdt);


	printf("\nTesting Square, Square Root, and Invert\n");
	for(i = 0; i < 6; i++)
	{
		a = ((double)rand()/(double)RAND_MAX)*exp(i);
		fa = DOUBLE_TO_FIX(a);

		//printf("%12f^2\t=%12f(double)\t%12f(fix)\terr: %.10f%%\n",
		printf("%12f,%12f,%12f,%.10f,",
			a, pow(a,2), FIX_TO_DOUBLE(fixsquare(fa)), 100.0*fixabs(pow(a,2)-(FIX_TO_DOUBLE(fixsquare(fa))))/pow(a,2));

		//printf("%10f^0.5\t=%12f(double)\t%12f(fix)\terr: %.10f%%\n",
		printf("%12f,%12f,%12f,%.10f,",
			a, sqrt(a), FIX_TO_DOUBLE(fixsqrt(fa)), 100.0*fixabs(sqrt(a)-(FIX_TO_DOUBLE(fixsqrt(fa))))/sqrt(a));

		//printf("1/%12f\t=%12f(double)\t%12f(fix)\terr: %.10f%%\n",
		printf("%12f,%12f,%12f,%.10f\n",
			a, 1/a, FIX_TO_DOUBLE(fixinv(fa)), 100.0*fixabs((1/a)-(FIX_TO_DOUBLE(fixinv(fa))))/fixabs(1/a));
	}

	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += dsqrt(dlist[j]+dlist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nSquare Root with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += fixsqrt(flist[j]+flist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	fdt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nSquare Root with fix: %12f s (%f), improvement %f\n", fdt, FIX_TO_DOUBLE(fa), dt/fdt);
	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += dinv(dlist[j]+dlist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nInvert with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += fixinv(flist[j]+flist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	fdt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nInvert with fix: %12f s (%f), improvement %f\n", fdt, FIX_TO_DOUBLE(fa), dt/fdt);

	printf("\nTesting Trigonometric Functions\n");
	for(i = -FIX_RAD_MAX+1; i <= FIX_RAD_MAX+1; i += FIX_RAD_MAX/16) //Offset avoids singularities
	{
		a = ((double)i)*2*M_PI/FIX_RAD_MAX;
		fa = fixdiv(fixmult(INT_TO_FIX(i),FIX_2PI),INT_TO_FIX(FIX_RAD_MAX));

		//printf("@%5d: sin(%10f)\t=%10f(double)\tsin(%10f)\t=%10f(fix)\terr: %.10f%%\n",
		printf("%5d,%10f,%10f,%10f,%10f,%.10f,",
			i, a, sin(a), FIX_TO_DOUBLE(fa), FIX_TO_DOUBLE(fixsin(fa)), 100.0*fixabs(sin(a)-(FIX_TO_DOUBLE(fixsin(fa))))/fixabs(sin(a)));

		//printf("@%5d: cos(%10f)\t=%10f(double)\tcos(%10f)\t=%10f(fix)\terr: %.10f%%\n",
		printf("%5d,%10f,%10f,%10f,%10f,%.10f\n",
			i, a, cos(a), FIX_TO_DOUBLE(fa), FIX_TO_DOUBLE(fixcos(fa)), 100.0*fixabs(cos(a)-(FIX_TO_DOUBLE(fixcos(fa))))/fixabs(cos(a)));

//		printf("@%5d: tan(%10f)\t=%10f(double)\ttan(%10f)\t=%10f(fix)\terr: %.10f%%\n",
		//printf("%5d,%10f,%10f,%10f,%10f,%.10f\n",
//			i, a, tan(a), FIX_TO_DOUBLE(fa), FIX_TO_DOUBLE(fixtan(fa)), 100.0*fixabs(tan(a)-(FIX_TO_DOUBLE(fixtan(fa))))/fixabs(tan(a)));
	}

	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += dsin(dlist[j]+dlist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nSine with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += fixsin(flist[j]+flist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nSine with fix: %12f s (%f), improvement %f\n", dt, FIX_TO_DOUBLE(fa), dt/fdt);
	printf("Timing Test\n");
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			a += dcos(dlist[j]+dlist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nCosine with double: %12f s (%f)\n", dt, a);
	gettimeofday(&timeStart, NULL);
	for(i = 0; i < TEST_LENGTH; i++)
		for(j = 0; j < TEST_LENGTH; j++) {
			fa += fixcos(flist[j]+flist[i]);
		}
	gettimeofday(&timeEnd, NULL);
	dt = (double)(timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0;
	printf("\nCosine with fix: %12f s (%f), improvement %f\n", dt, FIX_TO_DOUBLE(fa), dt/fdt);

	printf("\nMatrix Storage Order Tests\n");
	for(i = 0; i < MATRIX_SIZE; i++)
	{
		for(j = 0; j < MATRIX_SIZE; j++)
			printf("[%d][%d]@%lu=%f\t", i, j, (long unsigned int)&(mt[i][j]), FIX_TO_DOUBLE(mt[i][j]));
		printf("\n");
	}
	for(i = 0; i < MATRIX_SIZE*MATRIX_SIZE; i++)
		printf("[%d]@%lu=%f\n", i, (long unsigned int)(*mt+i), FIX_TO_DOUBLE(*(*mt+i)));

	printf("\nTesting Vector Operations\n");

	printf("\nSequential Vectors\n");
	for(j = 0; j < MATRIX_SIZE; j++)
		va[j] = INT_TO_FIX(j+1);
	PrintMatrix((fix*)va, 1, MATRIX_SIZE);
	fixvmake((fix*)vb, MATRIX_SIZE, DOUBLE_TO_FIX(0.1), DOUBLE_TO_FIX(0.2), DOUBLE_TO_FIX(0.3));
	PrintMatrix((fix*)vb, 1, MATRIX_SIZE);

	printf("\nVector Addition\n");
	fixvadd((fix*)vc, (fix*)va, (fix*)vb, MATRIX_SIZE);
	PrintMatrix((fix*)va, 1, MATRIX_SIZE);
	printf("+\n");
	PrintMatrix((fix*)vb, 1, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)vc, 1, MATRIX_SIZE);

	printf("\nVector Subtraction\n");
	fixvsub((fix*)vc, (fix*)va, (fix*)vb, MATRIX_SIZE);
	PrintMatrix((fix*)va, 1, MATRIX_SIZE);
	printf("+\n");
	PrintMatrix((fix*)vb, 1, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)vc, 1, MATRIX_SIZE);

	printf("\nVector Mutiplication (to matrix)\n");
	fixvmult((fix*)mc, (fix*)va, (fix*)(&mt[1]), MATRIX_SIZE);
	PrintMatrix((fix*)va, 1, MATRIX_SIZE);
	printf("*\n");
	PrintMatrix((fix*)(&mt[1]), 1, MATRIX_SIZE);
	printf("^T =\n");
	PrintMatrix((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nVector Dot Product\n");
	PrintMatrix((fix*)va, 1, MATRIX_SIZE);
	printf("*\n");
	PrintMatrix((fix*)vb, 1, MATRIX_SIZE);
	printf("= %f\n", FIX_TO_DOUBLE(fixvdot((fix*)va, (fix*)vb, MATRIX_SIZE)));

	printf("\nTesting Matrix Operations\n");

	printf("\nPositive-Definite Matrix\n");
	fixmmake((fix*)ma, MATRIX_SIZE, MATRIX_SIZE,
		INT_TO_FIX(2), INT_TO_FIX(-1), INT_TO_FIX(0),
		INT_TO_FIX(-1), INT_TO_FIX(2), INT_TO_FIX(-1),
		INT_TO_FIX(0), INT_TO_FIX(-1), INT_TO_FIX(2));
	PrintMatrix((fix*)ma, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nSequential Matrix\n");
	for(j = 0; j < MATRIX_SIZE*MATRIX_SIZE; j++)
		mb[j] = INT_TO_FIX(j+1);
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nZero Matrix\n");
	fixmzero((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);
	PrintMatrix((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nMatrix Addition\n");
	fixmadd((fix*)mc, (fix*)ma, (fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	PrintMatrix((fix*)ma, MATRIX_SIZE, MATRIX_SIZE);
	printf("+\n");
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nMatrix Subtraction\n");
	fixmsub((fix*)mc, (fix*)ma, (fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	PrintMatrix((fix*)ma, MATRIX_SIZE, MATRIX_SIZE);
	printf("-\n");
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nMatrix Multiplication by Square\n");
	fixmmultsq((fix*)mc, (fix*)ma, (fix*)mb, MATRIX_SIZE);
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("^T *\n");
	PrintMatrix((fix*)ma, MATRIX_SIZE, MATRIX_SIZE);
	printf("*\n");
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nMatrix Scalar Multiplication\n");
	PrintMatrix((fix*)mt, MATRIX_SIZE, MATRIX_SIZE);
	fixmscale((fix*)mt, DOUBLE_TO_FIX(MULTIPLIER), MATRIX_SIZE, MATRIX_SIZE);
	printf("* %f\n", MULTIPLIER);
	printf("=\n");
	PrintMatrix((fix*)mt, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nMatrix Vector Multiplication\n");
	fixvxform((fix*)vb, (fix*)mb, (fix*)va, MATRIX_SIZE, MATRIX_SIZE);
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("*\n");
	PrintMatrix((fix*)va, 1, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)vb, 1, MATRIX_SIZE);

	printf("\nMatrix Square Root\n");
	fixmchol((fix*)mb, (fix*)ma, MATRIX_SIZE);
	printf("CHOL\n");
	PrintMatrix((fix*)ma, MATRIX_SIZE, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nMatrix Transpose\n");
	fixmtrans((fix*)mc, (fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("^T\n");
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);

	printf("\nLU-Decomposed Matrix Squaring\n");
	fixmmult((fix*)ma, (fix*)mb, (fix*)mc, MATRIX_SIZE);
	PrintMatrix((fix*)mb, MATRIX_SIZE, MATRIX_SIZE);
	printf("*\n");
	PrintMatrix((fix*)mc, MATRIX_SIZE, MATRIX_SIZE);
	printf("=\n");
	PrintMatrix((fix*)ma, MATRIX_SIZE, MATRIX_SIZE);

	return 0;
}
