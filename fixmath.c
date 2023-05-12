//Based on "How to use Fixed Point (16.16) Math" by Night Stalker,
//http://netwinder.osuosl.org/pub/netwinder/docs/nw/fix1FAQ.html
//http://netwinder.osuosl.org/pub/netwinder/docs/nw/fix2FAQ.html
#include "../umath/fixmath.h"

#include <stdarg.h> //for va_* varargs

#include "../umath/fixtables.h"

//Inlining library functions: in C99, You have to put the definition (full body) in the header file, this will allow any file which includes the header file to be able to use the inline definition if the compiler chooses to do so.  You have to put the extern declaration (prototype) in the source file to tell the compiler to emit an extern version of the function in the library. This provides one place in your library for the non-inline version.

// Multiply
extern fix fixmult(fix n1, fix n2);
// Divide
extern fix fixdiv(fix num, fix den);
// Square
extern fix fixsquare(fix n);
// Invert
extern fix fixinv(fix n);
// Square Root
extern fix fixsqrt(fix n);
// Exponent
extern fix fixpow(fix n, int exp);
// Sine
extern fix fixsin(fix angle);
// Cosine
extern fix fixcos(fix angle);
// Tangent
extern fix fixtan(fix angle);

// Create Vector from Arguments
void fixvmake(fix *v, int length, ...)
{
	int i;
	va_list ap;
	va_start(ap, length);
	for(i = 0; i < length; i++)
		v[i] = va_arg(ap, fix);
	va_end(ap);
}

// Copy Vector
void fixvcopy(fix *v, fix *v1, int length)
{
	int i;
	for(i = 0; i < length; i++)
		v[i] =v1[i];
}

// Vector Zero
void fixvzero(fix *v, int length)
{
	int i;
	for(i = 0; i < length; i++)
		v[i] = FIX_ZERO;
}

// Vector One
void fixvone(fix *v, int length)
{
	int i;
	for(i = 0; i < length; i++)
		v[i] = FIX_ONE;
}

// Vector Scalar Addition
void fixvscalaradd(fix *v, fix *v1, fix s, int length)
{
	int i;
	for(i = 0; i < length; i++)
		v[i] = v1[i] + s;
}

// Vector Addition
void fixvadd(fix *v, fix *v1, fix *v2, int length)
{
	int i;
	for(i = 0; i < length; i++)
		v[i] = v1[i] + v2[i];
}

// Vector Subtraction
void fixvsub(fix *v, fix *v1, fix *v2, int length)
{
	int i;
	for(i = 0; i < length; i++)
		v[i] = v1[i] - v2[i];
}

// Vector Length
fix fixvlength(fix *v, int length)
{
	int i;
	fix l = FIX_ZERO;
	for(i = 0; i < length; i++)
		l += fixmult(v[i], v[i]);
	return fixsqrt(l);
}

// Vector Scalar Mutiplication
void fixvscale(fix *v, fix s, int length)
{
	int i;
	for(i = 0; i < length; i++)
		v[i] = fixmult(v[i], s);
}

// Vector Dot Product
fix fixvdot(fix *v1, fix *v2, int length)
{
	int i;
	fix vp = FIX_ZERO;
	for(i = 0; i < length; i++)
		vp += fixmult(v1[i], v2[i]);
	return vp;
}

// Vector Square Matrix Mutiplication (v1 * v2')
void fixvmult(fix *m, fix *v1, fix *v2, int length)
{
	int i, j;
	for(i = 0; i < length; i++)
		for(j = 0; j < length; j++)
			m[i*length+j] = fixmult(v1[i], v2[j]);
}

// Create Matrix from Arguments (row-major order)
void fixmmake(fix *m, int rows, int cols, ...)
{
	int i, j;
	va_list ap;
	va_start(ap, cols);
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = va_arg(ap, fix);
	va_end(ap);
}

// Create Diagonal Matrix from Arguments
void fixmmakediag(fix *m, int size, ...)
{
	int i;
	va_list ap;
	va_start(ap, size);
	for(i = 0; i < size; i++)
		m[i*size+i] = va_arg(ap, fix);
	va_end(ap);
}

// Copy Matrix
void fixmcopy(fix *m, fix *m1, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = m1[i*cols+j];
}

// Matrix Zero
void fixmzero(fix *m, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = FIX_ZERO;
}

// Matrix One
void fixmone(fix *m, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = FIX_ONE;
}

// Matrix Identity
void fixmeye(fix *m, int size)
{
	int i;
	for(i = 0; i < size; i++)
		m[i*size+i] = FIX_ONE;
}

// Matrix Diagonal
void fixmdiag(fix *m, fix *v1, int length)
{
	int i;
	for(i = 0; i < length; i++)
		m[i*length+i] = v1[i];
}

// Matrix Diagonal to Diagonal Matrix
void fixmdiagdiag(fix *m, fix *m1, int length)
{
	int i;
	for(i = 0; i < length; i++)
		m[i*length+i] = m1[i*length+i];
}

// Matrix Scalar Addition
void fixmscalaradd(fix *m, fix *m1, fix s, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = m1[i*cols+j] + s;
}

// Matrix Addition
void fixmadd(fix *m, fix *m1, fix *m2, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = m1[i*cols+j] + m2[i*cols+j];
}

// Matrix Subtraction
void fixmsub(fix *m, fix *m1, fix *m2, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = m1[i*cols+j] - m2[i*cols+j];
}

// Matrix Transposition
void fixmtrans(fix *m, fix *m1, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[j*cols+i] = m1[i*cols+j];
}

// Matrix Scalar Multiplication
void fixmscale(fix *m, fix s, int rows, int cols)
{
	int i, j;
	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			m[i*cols+j] = fixmult(m[i*cols+j],s);
}

// Matrix Vector Multiplication (transform)
void fixvxform(fix *v, fix *m1, fix *v2, int mrows, int mcols)
{
	int i, j;
	for(i = 0; i < mrows; i++)
	{
		v[i] = FIX_ZERO;
		for(j = 0; j < mcols; j++)
			v[i] = v[i] + fixmult(m1[i*mcols+j],v2[j]);
	}
}

// Matrix Multiplication (square matrices)
void fixmmult(fix *m, fix *m1, fix *m2, int size)
{
	int i, j, k;
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
		{
			m[i*size+j] = FIX_ZERO;
			for(k = 0; k < size; k++)
				m[i*size+j] = m[i*size+j] + fixmult(m1[i*size+k], m2[k*size+j]);
		}
}

// Matrix Multiplication by Square (M2 * M1 * M2' square matrices)
void fixmmultsq(fix *m, fix *m1, fix *m2, int size)
{
	int i, j, k;
	double x[size*size];
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
		{
			x[i*size+j] = 0;
			for(k = 0; k < size; k++)
				x[i*size+j] = x[i*size+j] + fixmult(m1[i*size+k], m2[j*size+k]);
		}
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
		{
			m[i*size+j] = 0;
			for(k = 0; k < size; k++)
				m[i*size+j] = m[i*size+j] + fixmult(m2[i*size+k], x[k*size+j]);
		}
}

// Matrix Square Root (square matrices)
void fixmchol(fix *msqrt, fix *m, int size)
{
	// Using Cholesky Decomposition
	int i, j, k;
	fix s;
	for (i = 0; i < size; i++) {
		for (j = 0; j < (i + 1); j++) {
			s = FIX_ZERO;
			for (k = 0; k < j; k++) {
				s += fixmult(msqrt[i*size+k], msqrt[j*size+k]);
			}
			msqrt[i*size+j] = ((i == j) ?
					fixsqrt(m[i*size+i] - s) :
					fixmult(fixdiv(FIX_ONE, msqrt[j*size+j]), (m[i*size+j] - s)));
		}
	}
}

// Vector Cross Product for 3D vectors
void fixv3cross(fix *v, fix *v1, fix *v2)
{
	v[0] = fixmult(v1[1], v2[2]) - fixmult(v1[2], v2[1]);
	v[1] = fixmult(v1[2], v2[0]) - fixmult(v1[0], v2[2]);
	v[2] = fixmult(v1[0], v2[1]) - fixmult(v1[1], v2[0]);
}

// Matrix Inverse of 2x2 matrix
int fixm2inv(fix *m, fix *m1)
{
	fix det = (fixmult(m1[0], m1[3]) - fixmult(m1[1], m1[2]));
	if(det == 0) return -1;
	m[0] = fixdiv(m1[3], det);
	m[1] = fixdiv(-m1[1], det);
	m[2] = fixdiv(-m1[2], det);
	m[3] = fixdiv(m1[0], det);
	return 0;
}

// Matrix Inverse of 3x3 matrix
int fixm3inv(fix *m, fix *m1)
{
	fix det = (
			fixmult(fixmult(m1[0],m1[4]),m1[8])
			- fixmult(fixmult(m1[0],m1[7]),m1[5])
			+ fixmult(fixmult(m1[3],m1[7]),m1[2])
			- fixmult(fixmult(m1[6],m1[4]),m1[2])
			+ fixmult(fixmult(m1[6],m1[1]),m1[5])
			- fixmult(fixmult(m1[3],m1[1]),m1[8]));
	if(det == 0) return -1;
	m[0] = fixdiv(fixmult(m1[4],m1[8]) - fixmult(m1[5],m1[7]), det);
	m[1] = fixdiv(fixmult(m1[2],m1[7]) - fixmult(m1[1],m1[8]), det);
	m[2] = fixdiv(fixmult(m1[1],m1[5]) - fixmult(m1[2],m1[4]), det);
	m[3] = fixdiv(fixmult(m1[5],m1[6]) - fixmult(m1[3],m1[8]), det);
	m[4] = fixdiv(fixmult(m1[0],m1[8]) - fixmult(m1[2],m1[6]), det);
	m[5] = fixdiv(fixmult(m1[2],m1[3]) - fixmult(m1[0],m1[5]), det);
	m[6] = fixdiv(fixmult(m1[3],m1[7]) - fixmult(m1[4],m1[6]), det);
	m[7] = fixdiv(fixmult(m1[1],m1[6]) - fixmult(m1[0],m1[7]), det);
	m[8] = fixdiv(fixmult(m1[0],m1[4]) - fixmult(m1[1],m1[3]), det);
	return 0;
}

// Matrix Inverse of 4x4 matrix
int fixm4inv(fix *m, fix *m1)
{
        fix det = (
          fixmult(fixmult(m1[0],m1[5]),fixmult(m1[10],m1[15])) + fixmult(fixmult(m1[0],m1[6]),fixmult(m1[11],m1[13])) + fixmult(fixmult(m1[0],m1[7]),fixmult(m1[9],m1[14]))
        - fixmult(fixmult(m1[0],m1[5]),fixmult(m1[11],m1[14])) + fixmult(fixmult(m1[0],m1[6]),fixmult(m1[9],m1[15])) + fixmult(fixmult(m1[0],m1[7]),fixmult(m1[10],m1[13]))
        + fixmult(fixmult(m1[1],m1[4]),fixmult(m1[11],m1[14])) + fixmult(fixmult(m1[1],m1[6]),fixmult(m1[8],m1[15])) + fixmult(fixmult(m1[1],m1[7]),fixmult(m1[10],m1[12]))
        - fixmult(fixmult(m1[1],m1[4]),fixmult(m1[10],m1[15])) + fixmult(fixmult(m1[1],m1[6]),fixmult(m1[11],m1[12])) + fixmult(fixmult(m1[1],m1[7]),fixmult(m1[8],m1[14]))
        + fixmult(fixmult(m1[2],m1[4]),fixmult(m1[9],m1[15])) + fixmult(fixmult(m1[2],m1[5]),fixmult(m1[11],m1[12])) + fixmult(fixmult(m1[2],m1[7]),fixmult(m1[8],m1[13]))
        - fixmult(fixmult(m1[2],m1[4]),fixmult(m1[11],m1[13])) + fixmult(fixmult(m1[2],m1[5]),fixmult(m1[8],m1[15])) + fixmult(fixmult(m1[2],m1[7]),fixmult(m1[9],m1[12]))
        + fixmult(fixmult(m1[3],m1[4]),fixmult(m1[10],m1[13])) + fixmult(fixmult(m1[3],m1[5]),fixmult(m1[8],m1[14])) + fixmult(fixmult(m1[3],m1[6]),fixmult(m1[9],m1[12]))
        - fixmult(fixmult(m1[3],m1[4]),fixmult(m1[9],m1[14])) + fixmult(fixmult(m1[3],m1[5]),fixmult(m1[10],m1[12])) + fixmult(fixmult(m1[3],m1[6]),fixmult(m1[8],m1[13])));
        if(det == 0) return -1;
        m[0] = fixdiv(fixmult(fixmult(m1[5],m1[10]),m1[15]) + fixmult(fixmult(m1[6],m1[11]),m1[13]) + fixmult(fixmult(m1[7],m1[9]),m1[14]) - fixmult(fixmult(m1[5],m1[11]),m1[14]) - fixmult(fixmult(m1[6],m1[9]),m1[15]) - fixmult(fixmult(m1[7],m1[10]),m1[13]), det);
        m[1] = fixdiv(fixmult(fixmult(m1[1],m1[11]),m1[14]) + fixmult(fixmult(m1[2],m1[9]),m1[15]) + fixmult(fixmult(m1[3],m1[10]),m1[13]) - fixmult(fixmult(m1[1],m1[12]),m1[15]) - fixmult(fixmult(m1[2],m1[11]),m1[13]) - fixmult(fixmult(m1[3],m1[9]),m1[14]), det);
        m[2] = fixdiv(fixmult(fixmult(m1[1],m1[6]),m1[15]) + fixmult(fixmult(m1[2],m1[7]),m1[13]) + fixmult(fixmult(m1[3],m1[5]),m1[14]) - fixmult(fixmult(m1[5],m1[7]),m1[14]) - fixmult(fixmult(m1[2],m1[5]),m1[15]) - fixmult(fixmult(m1[3],m1[6]),m1[13]), det);
        m[3] = fixdiv(fixmult(fixmult(m1[1],m1[7]),m1[10]) + fixmult(fixmult(m1[2],m1[5]),m1[11]) + fixmult(fixmult(m1[3],m1[6]),m1[9]) - fixmult(fixmult(m1[5],m1[6]),m1[11]) - fixmult(fixmult(m1[2],m1[7]),m1[9]) - fixmult(fixmult(m1[3],m1[5]),m1[10]), det);
        m[4] = fixdiv(fixmult(fixmult(m1[4],m1[11]),m1[14]) + fixmult(fixmult(m1[6],m1[8]),m1[15]) + fixmult(fixmult(m1[7],m1[10]),m1[12]) - fixmult(fixmult(m1[4],m1[10]),m1[15]) - fixmult(fixmult(m1[6],m1[11]),m1[12]) - fixmult(fixmult(m1[7],m1[8]),m1[14]), det);
        m[5] = fixdiv(fixmult(fixmult(m1[0],m1[10]),m1[15]) + fixmult(fixmult(m1[2],m1[11]),m1[12]) + fixmult(fixmult(m1[3],m1[8]),m1[14]) - fixmult(fixmult(m1[0],m1[11]),m1[14]) - fixmult(fixmult(m1[2],m1[8]),m1[15]) - fixmult(fixmult(m1[3],m1[10]),m1[12]), det);
        m[6] = fixdiv(fixmult(fixmult(m1[0],m1[7]),m1[14]) + fixmult(fixmult(m1[2],m1[4]),m1[15]) + fixmult(fixmult(m1[3],m1[6]),m1[12]) - fixmult(fixmult(m1[0],m1[6]),m1[15]) - fixmult(fixmult(m1[2],m1[7]),m1[12]) - fixmult(fixmult(m1[3],m1[4]),m1[14]), det);
        m[7] = fixdiv(fixmult(fixmult(m1[0],m1[6]),m1[11]) + fixmult(fixmult(m1[2],m1[7]),m1[8]) + fixmult(fixmult(m1[3],m1[4]),m1[10]) - fixmult(fixmult(m1[0],m1[7]),m1[10]) - fixmult(fixmult(m1[2],m1[4]),m1[11]) - fixmult(fixmult(m1[3],m1[6]),m1[8]), det);
        m[8] = fixdiv(fixmult(fixmult(m1[4],m1[9]),m1[15]) + fixmult(fixmult(m1[5],m1[11]),m1[12]) + fixmult(fixmult(m1[7],m1[8]),m1[13]) - fixmult(fixmult(m1[4],m1[11]),m1[13]) - fixmult(fixmult(m1[5],m1[8]),m1[15]) - fixmult(fixmult(m1[7],m1[9]),m1[12]), det);
        m[9] = fixdiv(fixmult(fixmult(m1[0],m1[11]),m1[13]) + fixmult(fixmult(m1[1],m1[8]),m1[15]) + fixmult(fixmult(m1[3],m1[9]),m1[12]) - fixmult(fixmult(m1[0],m1[9]),m1[15]) - fixmult(fixmult(m1[1],m1[11]),m1[12]) - fixmult(fixmult(m1[3],m1[8]),m1[13]), det);
        m[10] =fixdiv(fixmult(fixmult(m1[0],m1[5]),m1[15]) + fixmult(fixmult(m1[1],m1[7]),m1[12]) + fixmult(fixmult(m1[3],m1[4]),m1[13]) - fixmult(fixmult(m1[0],m1[7]),m1[13]) - fixmult(fixmult(m1[1],m1[4]),m1[15]) - fixmult(fixmult(m1[3],m1[5]),m1[12]), det);
        m[11] =fixdiv(fixmult(fixmult(m1[0],m1[7]),m1[9]) + fixmult(fixmult(m1[1],m1[4]),m1[11]) + fixmult(fixmult(m1[3],m1[5]),m1[8]) - fixmult(fixmult(m1[0],m1[5]),m1[11]) - fixmult(fixmult(m1[1],m1[7]),m1[8]) - fixmult(fixmult(m1[3],m1[4]),m1[9]), det);
        m[12] =fixdiv(fixmult(fixmult(m1[4],m1[10]),m1[13]) + fixmult(fixmult(m1[5],m1[8]),m1[14]) + fixmult(fixmult(m1[6],m1[9]),m1[12]) - fixmult(fixmult(m1[4],m1[9]),m1[14]) - fixmult(fixmult(m1[5],m1[10]),m1[12]) - fixmult(fixmult(m1[6],m1[8]),m1[13]), det);
        m[13] =fixdiv(fixmult(fixmult(m1[0],m1[9]),m1[14]) + fixmult(fixmult(m1[1],m1[10]),m1[12]) + fixmult(fixmult(m1[2],m1[8]),m1[13]) - fixmult(fixmult(m1[0],m1[10]),m1[13]) - fixmult(fixmult(m1[1],m1[8]),m1[14]) - fixmult(fixmult(m1[2],m1[9]),m1[12]), det);
        m[14] =fixdiv(fixmult(fixmult(m1[0],m1[6]),m1[13]) + fixmult(fixmult(m1[1],m1[4]),m1[14]) + fixmult(fixmult(m1[2],m1[5]),m1[12]) - fixmult(fixmult(m1[0],m1[5]),m1[14]) - fixmult(fixmult(m1[1],m1[6]),m1[12]) - fixmult(fixmult(m1[2],m1[4]),m1[13]), det);
        m[15] =fixdiv(fixmult(fixmult(m1[0],m1[5]),m1[10]) + fixmult(fixmult(m1[1],m1[6]),m1[8]) + fixmult(fixmult(m1[2],m1[4]),m1[9]) - fixmult(fixmult(m1[0],m1[6]),m1[9]) - fixmult(fixmult(m1[1],m1[4]),m1[10]) - fixmult(fixmult(m1[2],m1[5]),m1[8]), det);
        return 0;
}

// EigenValues of 2X2 Square Matrix
void dm2eigenval(fix *evalue, fix *m)
{
	//TODO: Should these be reversed in order?  It is this way to match with MATLAB
	evalue[3] = fixdiv(m[0]+m[3],INT_TO_FIX(2)) + fixsqrt(fixdiv(fixmult(m[0]+m[3],m[0]+m[3]),INT_TO_FIX(4)) + fixmult(m[1],m[2])-fixmult(m[0],m[3]));
	evalue[0] = fixdiv(m[0]+m[3],INT_TO_FIX(2)) - fixsqrt(fixdiv(fixmult(m[0]+m[3],m[0]+m[3]),INT_TO_FIX(4)) + fixmult(m[1],m[2])-fixmult(m[0],m[3]));
	evalue[1] = FIX_ZERO;
	evalue[2] = FIX_ZERO;
	// Same values, but different calculation
	//evalue[1] = (m[0]+m[3])/2.0 + sqrt(4.0*m[1]*m[2] + (m[0]-m[3])*(m[0]-m[3]))/2.0;
	//evalue[2] = (m[0]+m[3])/2.0 - sqrt(4.0*m[1]*m[2] + (m[0]-m[3])*(m[0]-m[3]))/2.0;
}

// EigenVectors from EigenValues of 2X2 Square Matrix
void dm2eigenvect(fix *evect, fix *m, fix *evalue)
{
	//TODO: Should these be negated?  It is this way to match with MATLAB
	evect[0] = fixdiv(-FIX_ONE,fixsqrt(FIX_ONE+fixpow(fixdiv(evalue[0]-m[0],m[1]),2)));
	evect[2] = -fixdiv(fixdiv(evalue[0]-m[0],m[1]),fixsqrt(1+fixpow(fixdiv(evalue[0]-m[0],m[1]),2)));
	evect[1] = fixdiv(-FIX_ONE,fixsqrt(FIX_ONE+fixpow(fixdiv(evalue[3]-m[0],m[1]),2)));
	evect[3] = -fixdiv(fixdiv(evalue[3]-m[0],m[1]),fixsqrt(1+fixpow(fixdiv(evalue[3]-m[0],m[1]),2)));
}

// 3D Vector Functions
void fixvect_make(fvector v, fix x, fix y, fix z)
{
	*v++ = x;
	*v++ = y;
	*v   = z;
}

void fixvect_copy(fvector sv, fvector dv)
{
	*dv++ = *sv++;
	*dv++ = *sv++;
	*dv   = *sv;
}

void fixvect_add(fvector v1, fvector v2, fvector dv)
{
	*dv++ = v1[VECT_X] + v2[VECT_X];
	*dv++ = v1[VECT_Y] + v2[VECT_Y];
	*dv   = v1[VECT_Z] + v2[VECT_Z];
}

void fixvect_sub(fvector v1, fvector v2, fvector dv)
{
	*dv++ = v1[VECT_X] - v2[VECT_X];
	*dv++ = v1[VECT_Y] - v2[VECT_Y];
	*dv   = v1[VECT_Z] - v2[VECT_Z];
}

fix fixvect_len(fvector v)
{
	// We use the High Precision Square Root here.
	return fixsqrt(fixsquare(*v) + fixsquare(v[1]) + fixsquare(v[2]));
}

fix fixvect_len2(fvector v)
{
	return (fixsquare(*v) + fixsquare(v[1]) + fixsquare(v[2]));
}

fix fixvect_dist(fvector v1, fvector v2)
{
	return fixsqrt(	fixsquare((*v1) - (*v2)) +
			fixsquare(v1[1] - v2[1]) +
			fixsquare(v1[2] - v2[2]));
}

fix fixvect_dist2(fvector v1, fvector v2)
{
	return (fixsquare((*v1) - (*v2)) + fixsquare(v1[1] - v2[1]) +
			fixsquare(v1[2] - v2[2]));
}

void fixvect_scale(fvector v, fix factor)
{
	*v = fixmult(*v, factor);
	v++;
	*v = fixmult(*v, factor);
	v++;
	*v = fixmult(*v, factor);
}

void fixvect_normalize(fvector v)
{
	fix factor = fixinv(fixvect_len(v));

	*v = fixmult(*v, factor);
	v++;
	*v = fixmult(*v, factor);
	v++;
	*v = fixmult(*v, factor);
}

fix fixvect_dot(fvector v1, fvector v2)
{
	return (fixmult(v1[VECT_X], v2[VECT_X]) +
			fixmult(v1[VECT_Y], v2[VECT_Y]) +
			fixmult(v1[VECT_Z], v2[VECT_Z]));
}

void fixvect_cross(fvector dv, fvector v1, fvector v2)
{
	*dv++ = fixmult(v1[VECT_Y], v2[VECT_Z]) - fixmult(v1[VECT_Z], v2[VECT_Y]);
	*dv++ = fixmult(v1[VECT_Z], v2[VECT_X]) - fixmult(v1[VECT_X], v2[VECT_Z]);
	*dv   = fixmult(v1[VECT_X], v2[VECT_Y]) - fixmult(v1[VECT_Y], v2[VECT_X]);
}

fix fixvect_crossz(fvector v1, fvector v2)
{
	return fixmult(*v1, v2[VECT_Y]) - fixmult(v1[VECT_Y], *v2);
}

// 3D Matrix Functions
void fixmatrix_make_zero(fmatrix m)
{
	m[0] = FIX_ZERO;
	m[1] = FIX_ZERO;
	m[2] = FIX_ZERO;
	m[3] = FIX_ZERO;
	m[4] = FIX_ZERO;
	m[5] = FIX_ZERO;
	m[6] = FIX_ZERO;
	m[7] = FIX_ZERO;
	m[8] = FIX_ZERO;
	m[9] = FIX_ZERO;
	m[10] = FIX_ZERO;
	m[11] = FIX_ZERO;
}

void fixmatrix_make_one(fmatrix m)
{
	m[0] = FIX_ONE;
	m[1] = FIX_ONE;
	m[2] = FIX_ONE;
	m[3] = FIX_ONE;
	m[4] = FIX_ONE;
	m[5] = FIX_ONE;
	m[6] = FIX_ONE;
	m[7] = FIX_ONE;
	m[8] = FIX_ONE;
	m[9] = FIX_ONE;
	m[10] = FIX_ONE;
	m[11] = FIX_ONE;
}

void fixmatrix_make_eye(fmatrix m)
{
	m[0] = FIX_ONE;
	m[1] = FIX_ZERO;
	m[2] = FIX_ZERO;
	m[3] = FIX_ZERO;
	m[4] = FIX_ZERO;
	m[5] = FIX_ONE;
	m[6] = FIX_ZERO;
	m[7] = FIX_ZERO;
	m[8] = FIX_ZERO;
	m[9] = FIX_ZERO;
	m[10] = FIX_ONE;
	m[11] = FIX_ZERO;
}

void fixmatrix_make_scale(fmatrix m, fvector sv)
{
	*m++ = *sv++;  *m++ =0;  *m++ =  0;   *m++ =  0;
	*m++ =0;  *m++ = *sv++;  *m++ =  0;   *m++ =  0;
	*m++ =0;  *m++ =0;  *m++ = *sv;  *m   =  0;
}

void fixmatrix_make_trans(fmatrix m, fvector tv)
{
	*m++ = FIX_ONE;  *m++ =   0;  *m++ =   0;  *m++ = *tv++;
	*m++ =   0;  *m++ = FIX_ONE;  *m++ =   0;  *m++ = *tv++;
	*m++ =   0;  *m++ =   0;  *m++ = FIX_ONE;  *m   = *tv;
}

void fixmatrix_make_xrot(fmatrix m, fixrad theta)
{
	fix trig_sin = fixsin(theta), trig_cos = fixcos(theta);
	*m++ = FIX_ONE;  *m++ =	0;  *m++ =	0;  *m++ = 0;
	*m++ =   0;  *m++ =  trig_cos;  *m++ =  trig_sin;  *m++ = 0;
	*m++ =   0;  *m++ = -trig_sin;  *m++ = -trig_cos;  *m   = 0;
}

void fixmatrix_make_yrot(fmatrix m, fixrad theta)
{
	fix trig_sin = fixsin(theta), trig_cos = fixcos(theta);
	*m++ =  trig_cos;  *m++ =   0;  *m++ = -trig_sin;  *m++ = 0;
	*m++ =	0;  *m++ = FIX_ONE;  *m++ =	0;  *m++ = 0;
	*m++ = -trig_sin;  *m++ =   0;  *m++ =  trig_cos;  *m   = 0;
}

void fixmatrix_make_zrot(fmatrix m, fixrad theta)
{
	fix trig_sin = fixsin(theta), trig_cos = fixcos(theta);
	*m++ =  trig_cos;  *m++ =  trig_sin;  *m++ =   0;  *m++ = 0;
	*m++ = -trig_sin;  *m++ =  trig_cos;  *m++ =   0;  *m++ = 0;
	*m++ =	0;  *m++ =	0;  *m++ = FIX_ONE;  *m   = 0;
}

void fixmatrix_copy(fmatrix dm, fmatrix sm)
{
	int i;
	for(i = 0; i < 12; i++)
		dm[i] = sm[i];
}

void fixmatrix_concatinate(fmatrix dm, fmatrix m1, fmatrix m2)
{
	int i = 3, temp = 0;
	while (i--) {
		*dm++ = fixmult(m1[temp  ], m2[0]) +
				fixmult(m1[temp+1], m2[4]) +
				fixmult(m1[temp+2], m2[8]);

		*dm++ = fixmult(m1[temp  ], m2[1]) +
				fixmult(m1[temp+1], m2[5]) +
				fixmult(m1[temp+2], m2[9]);

		*dm++ = fixmult(m1[temp  ], m2[2]) +
				fixmult(m1[temp+1], m2[6]) +
				fixmult(m1[temp+2], m2[10]);

		*dm++ = fixmult(m1[temp  ], m2[3]) +
				fixmult(m1[temp+1], m2[7]) +
				fixmult(m1[temp+2], m2[11]) + m1[temp+3];

		temp += 4;
	}
}

void fixmatrix_trans(fmatrix m, fvector dv)
{
	m[ 3] += *dv++;
	m[ 7] += *dv++;
	m[11] += *dv;
}

void fixmatrix_append_xrot(fmatrix m, fixrad theta)
{
	fix trig_sin = fixsin(theta), trig_cos = fixcos(theta);
	fix temp0, temp1, temp2;
	temp0 = fixmult(trig_cos, m[4]) + fixmult(-trig_sin, m[8]);
	temp1 = fixmult(trig_cos, m[5]) + fixmult(-trig_sin, m[9]);
	temp2 = fixmult(trig_cos, m[6]) + fixmult(-trig_sin, m[10]);

	m[8]  = fixmult(trig_sin, m[4]) + fixmult(trig_cos, m[8]);
	m[9]  = fixmult(trig_sin, m[5]) + fixmult(trig_cos, m[9]);
	m[10] = fixmult(trig_sin, m[6]) + fixmult(trig_cos, m[10]);

	m[4]  = temp0;
	m[5]  = temp1;
	m[6]  = temp2;
}

void fixmatrix_append_yrot(fmatrix m, fixrad theta)
{
	fix trig_sin = fixsin(theta), trig_cos = fixcos(theta);
	fix temp0, temp1, temp2;
	temp0 = fixmult(trig_cos, m[0]) + fixmult(trig_sin, m[8]);
	temp1 = fixmult(trig_cos, m[1]) + fixmult(trig_sin, m[9]);
	temp2 = fixmult(trig_cos, m[2]) + fixmult(trig_sin, m[10]);

	m[8]  = fixmult(-trig_sin, m[0]) + fixmult(trig_cos, m[8]);
	m[9]  = fixmult(-trig_sin, m[1]) + fixmult(trig_cos, m[9]);
	m[10] = fixmult(-trig_sin, m[2]) + fixmult(trig_cos, m[10]);

	m[0]  = temp0;
	m[1]  = temp1;
	m[2]  = temp2;
}

void fixmatrix_append_zrot(fmatrix m, fixrad theta)
{
	fix trig_sin = fixsin(theta), trig_cos = fixcos(theta);
	fix temp0, temp1, temp2;
	temp0 = fixmult(trig_cos, m[0]) + fixmult(-trig_sin, m[4]);
	temp1 = fixmult(trig_cos, m[1]) + fixmult(-trig_sin, m[5]);
	temp2 = fixmult(trig_cos, m[2]) + fixmult(-trig_sin, m[6]);

	m[4] = fixmult(trig_sin, m[0]) + fixmult(trig_cos, m[4]);
	m[5] = fixmult(trig_sin, m[1]) + fixmult(trig_cos, m[5]);
	m[6] = fixmult(trig_sin, m[2]) + fixmult(trig_cos, m[6]);

	m[0] = temp0;
	m[1] = temp1;
	m[2] = temp2;
}

void fixmatrix_xform(fvector dv, fvector sv, fmatrix m)
{
	fix sv0, sv1, sv2;
	sv0 = *sv++;
	sv1 = *sv++;
	sv2 = *sv;

	*dv = fixmult(*m++, sv0);
	*dv += fixmult(*m++, sv1);
	*dv += fixmult(*m++, sv2);
	*dv++ += *m++;

	*dv = fixmult(*m++, sv0);
	*dv += fixmult(*m++, sv1);
	*dv += fixmult(*m++, sv2);
	*dv++ += *m++;

	*dv = fixmult(*m++, sv0);
	*dv += fixmult(*m++, sv1);
	*dv += fixmult(*m++, sv2);
	*dv += *m;
}
