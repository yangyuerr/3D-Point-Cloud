

#ifndef _SVD_H_
#define _SVD_H_

#ifdef __cplusplus
extern "C" {
#endif

#define M 3
#define N 3

#include <stdio.h>
#include <stdlib.h>


#include "math.h"

#define NR_END 1
#define FREE_ARG char*

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define FMAX(a,b) ((a)>(b)? (a) : (b))
#define SQR(x) ((x)*(x))
#define IMIN(a,b) ((a) < (b) ?(a) : (b))

	/* Numerical Recipes standard error handler */
	void nrerror(char error_text[]);

	/* allocate a double vector with subscript range v[nl..nh] */
	double *fvector(long nl, long nh);

	/* free a double vector allocated with vector() */
	void free_fvector(double *v, long nl, long nh);

	/* svd operation api */
	double pythag(double a, double b);
	void svdcmp(double **a, int m, int n, double w[], double **v);

	void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);

	/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
	double **fmatrix(long nrl, long nrh, long ncl, long nch);

	/* free a double matrix allocated by matrix() */
	void free_fmatrix(double **m, long nrl, long nrh, long ncl, long nch);

	double get_norm(double *x, int n);
	double normalize(double *x, int n);
	double product(double*a, double *b, int n);
	void orth(double *a, double *b, int n);

/*
A: 输入待分解矩阵
K: 输入，取前K大奇异值及奇异向量
U[0],U[1],...,U[K-1]: 前K大奇异值对应的左奇异向量
S[0],S[1],...,S[K-1]: 前K大奇异值 S[0]>=S[1]>=...>=S[K-1]
V[0],V[1],...,V[K-1]: 前K大奇异值对应的右奇异向量
*/
bool svd(double A[][N], int K, double U[][N], double S[N], double V[][N]);

#ifdef __cplusplus
}
#endif

#endif