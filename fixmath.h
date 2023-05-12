//Based on "How to use Fixed Point (16.16) Math" by Night Stalker,
//http://netwinder.osuosl.org/pub/netwinder/docs/nw/fix1FAQ.html
//http://netwinder.osuosl.org/pub/netwinder/docs/nw/fix2FAQ.html
#pragma once

#include <inttypes.h>
#include <limits.h>

//Fixed Point type definitions
//Maximum value for 16.16: 65535+65535/65536 = 65535.999984741211
typedef int32_t fix;		// 16.16 FixedPoint
typedef uint32_t ufix;		// Unsigned type for square root
typedef int64_t longfix;	// Longer type for operations
typedef uint16_t fixrad;	// Integer angle (0..1023)
typedef fix fvector[3];		// Vector
typedef fix fmatrix[12];	// Matrix[3][4]
//Number of bits to the right of the decimal point
#define FIX_FRAC_BITS		16
#define FIX_SCALE		(1<<FIX_FRAC_BITS)
#define FIX_STEP		(1.0/(float)FIX_SCALE)
//Maximum number that can be represented in fixed point notation
#define FIX_MAX_WHOLE		((1<<(sizeof(fix)*8-FIX_FRAC_BITS))-1)
#define FIX_MAX_FRACTION	(((float)FIX_SCALE-1)/(float)FIX_SCALE)
#define FIX_MAX			((double)FIX_MAX_WHOLE + (double)FIX_MAX_FRACTION)
//Number of distinct angles in a revolution
#define FIX_RAD_MAX		1024
#define FIX_RAD_SCALE		162.974

//Constants
#define FIX_ZERO		0L
#define FIX_ONE			FIX_SCALE
#define FIX_PI			205887L
#define FIX_2PI			411775L
#define FIX_E			178144L
#define FIX_ROOT2		74804L
#define FIX_ROOT3		113512L
#define FIX_GOLDEN		106039L

//Macros for 3D Vector access
#define VECT_X	0
#define VECT_Y	1
#define VECT_Z	2

//Conversion Macros
#define INT_TO_FIX(x)		(fix)((x) << FIX_FRAC_BITS)
#define DOUBLE_TO_FIX(x)	(fix)((longfix)(x * (longfix)(FIX_SCALE + 0.5)))
#define FIX_TO_INT(x)		(int)((x) >> FIX_FRAC_BITS)
#define FIX_TO_DOUBLE(x)	(((double)x) / (double)FIX_SCALE)
#define FIX_TO_FLOAT(x)		(((float)x) / (float)FIX_SCALE)
#define ROUND_FIX_TO_INT(x)	(int)(((x) + 1<<(FIX_FRAC_BITS-1)) >> FIX_FRAC_BITS)
#define ROUND_FIX_TO_RAD(x)	(fixrad)(((fixmult(x,DOUBLE_TO_FIX(FIX_RAD_SCALE))) + (1<<(FIX_FRAC_BITS-1))) >> FIX_FRAC_BITS)

//Utility Macros
#define VSIZE(array) (sizeof(array) / sizeof((array)[0]))

// Math Operations
// Absolute Value
inline fix fixabs(fix n) { return ((n^(n>>(sizeof(fix)-1))) - (n>>(sizeof(fix)-1))); }

// Multiplication
inline fix fixmult(fix n1, fix n2) { return (fix)(((longfix)n1*(longfix)n2) >> FIX_FRAC_BITS); }

// Division
inline fix fixdiv(fix num, fix den) { return (fix)((((longfix)num << (FIX_FRAC_BITS)) / (longfix)den)); }

// This is faster than using mult for squares
inline fix fixsquare(fix n) { return (fix)(((longfix)n*(longfix)n) >> FIX_FRAC_BITS); }

// This is faster than using div
inline fix fixinv(fix n) { return (fix)((((longfix)FIX_ONE << (FIX_FRAC_BITS)) / (longfix)n)); }

// Square Root
inline fix fixsqrt(fix n)
{
	ufix root, remHi, remLo, testDiv, count;
	root = 0;	/* Clear root */
	remHi = 0;	/* Clear high part of partial remainder */
	remLo = n;	/* Get argument into low part of partial remainder */
	count = (15 + (FIX_FRAC_BITS >> 1)); /* Load loop counter, Must be even! */

	do {
		remHi = (remHi << 2) | (remLo >> 30); remLo <<= 2;	/* get 2 bits of arg */
		root <<= 1;	/* Get ready for the next bit in the root */
		testDiv = (root << 1) + 1;	/* Test radical */
		if (remHi >= testDiv) {
			remHi -= testDiv;
			root += 1;
		}
	} while (count-- != 0);

	return(root);
}

// Exponent
inline fix fixpow(fix n, int exp)
{
	while(exp-- > 1)
		n = fixmult(n, n);
	return n;
}

// Trigonometric Functions
extern fix SinTable[], CosTable[];
inline fix fixsin(fix angle)
{
	fixrad index = ROUND_FIX_TO_RAD(angle);
	index &= (FIX_RAD_MAX - 1);
	return SinTable[index];
}

inline fix fixcos(fix angle)
{
	fixrad index = ROUND_FIX_TO_RAD(angle);
	index &= (FIX_RAD_MAX - 1);
	return CosTable[index];
}

inline fix fixtan(fix angle)
{
	fixrad index = ROUND_FIX_TO_RAD(angle);
	index &= (FIX_RAD_MAX - 1);
	return (fixdiv(SinTable[index] << 16, CosTable[index]) >> 16);
}

// Vector Operations
// Create Vector from Arguments
void fixvmake(fix *v, int length, ...);
// Copy Vector
void fixvcopy(fix *v, fix *v1, int length);
// Vector Zero
void fixvzero(fix *v, int length);
// Vector One
void fixvone(fix *v, int length);
// Vector Scalar Addition
void fixvscalaradd(fix *v, fix *v1, fix s, int length);
// Vector Addition
void fixvadd(fix *v, fix *v1, fix *v2, int length);
// Vector Subtraction
void fixvsub(fix *v, fix *v1, fix *v2, int length);
// Vector Length
fix fixvlength(fix *v, int length);
// Vector Scalar Multiplication
void fixvscale(fix *v, fix s, int length);
// Vector Dot Product
fix fixvdot(fix *v1, fix *v2, int length);
// Vector Square Matrix Multiplication (v1 * v2')
void fixvmult(fix *m, fix *v1, fix *v2, int length);

// Matrix Operations
// Create Matrix from Arguments (row-major order)
void fixmmake(fix *m, int rows, int cols, ...);
// Create Diagonal Matrix from Arguments
void fixmmakediag(fix *m, int size, ...);
// Copy Matrix
void fixmcopy(fix *m, fix *m1, int rows, int cols);
// Matrix Zero
void fixmzero(fix *m, int rows, int cols);
// Matrix One
void fixmone(fix *m, int rows, int cols);
// Matrix Identity
void fixmeye(fix *m, int size);
// Matrix Diagonal
void fixmdiag(fix *m, fix *v1, int length);
// Matrix Diagonal to Diagonal Matrix
void fixmdiagdiag(fix *m, fix *m1, int length);
// Matrix Scalar Addition
void fixmscalaradd(fix *m, fix *m1, fix s, int rows, int cols);
// Matrix Addition
void fixmadd(fix *m, fix *m1, fix *m2, int rows, int cols);
// Matrix Subtraction
void fixmsub(fix *m, fix *m1, fix *m2, int rows, int cols);
// Matrix Transposition
void fixmtrans(fix *m, fix *m1, int rows, int cols);
// Matrix Scalar Multiplication
void fixmscale(fix *m, fix s, int rows, int cols);

// Algebraic Operations
// Matrix Vector Multiplication (transform)
void fixvxform(fix *v, fix *m1, fix *v2, int mrows, int mcols);
// Matrix Multiplication (square matrices)
void fixmmult(fix *m, fix *m1, fix *m2, int size);
// Matrix Multiplication by Square (M2 * M1 * M2' square matrices)
void fixmmultsq(fix *m, fix *m1, fix *m2, int size);
// Matrix Square Root (square matrices)
void fixmchol(fix *m, fix *m1, int size);
// Vector Cross Product for 3D vectors
void fixv3cross(fix *v, fix *v1, fix *v2);
// Matrix Inverse of 2x2 matrix
int fixm2inv(fix *m, fix *m1);
// Matrix Inverse of 3x3 matrix
int fixm3inv(fix *m, fix *m1);
// Matrix Inverse of 4x4 matrix
int fixm4inv(fix *m, fix *m1);
// EigenValues of 2X2 Square Matrix
void fixm2eigenval(fix *evalue, fix *m);
// EigenVectors from EigenValues of 2X2 Square Matrix
void fixm2eigenvect(fix *evect, fix *m, fix *evalue);

// 3D Vector Functions
void fixvect_make(fvector v, fix x, fix y, fix z);
void fixvect_copy(fvector sv, fvector dv);
void fixvect_add(fvector v1, fvector v2, fvector dv);
void fixvect_sub(fvector v1, fvector v2, fvector dv);
fix fixvect_len(fvector v);
fix fixvect_len2(fvector v);
fix fixvect_dist(fvector v1, fvector v2);
fix fixvect_dist2(fvector v1, fvector v2);
void fixvect_scale(fvector v, fix factor);
void fixvect_normalize(fvector v);
fix fixvect_dot(fvector v1, fvector v2);
void fixvect_cross(fvector dv, fvector v1, fvector v2);
fix fixvect_crossz(fvector v1, fvector v2);

// 3D Matrix Functions
void fixmatrix_make_zero(fmatrix m);
void fixmatrix_make_one(fmatrix m);
void fixmatrix_make_eye(fmatrix m);
void fixmatrix_make_scale(fmatrix m, fvector sv);
void fixmatrix_make_trans(fmatrix m, fvector tv);
void fixmatrix_make_xrot(fmatrix m, fixrad theta);
void fixmatrix_make_yrot(fmatrix m, fixrad theta);
void fixmatrix_make_zrot(fmatrix m, fixrad theta);
void fixmatrix_copy(fmatrix dm, fmatrix sm);
void fixmatrix_concatinate(fmatrix dm, fmatrix m1, fmatrix m2);
void fixmatrix_trans(fmatrix m, fvector dv);
void fixmatrix_append_xrot(fmatrix m, fixrad theta);
void fixmatrix_append_yrot(fmatrix m, fixrad theta);
void fixmatrix_append_zrot(fmatrix m, fixrad theta);
void fixmatrix_xform(fvector dv, fvector sv, fmatrix m);
