/*
bool svd(vector<vector<double> > A, int K, vector<vector<double> > &U, vector<double> &S, vector<vector<double> > &V);
A: 输入待分解矩阵
K: 输入，取前K大奇异值及奇异向量
U[0],U[1],...,U[K-1]: 前K大奇异值对应的左奇异向量
S[0],S[1],...,S[K-1]: 前K大奇异值 S[0]>=S[1]>=...>=S[K-1]
V[0],V[1],...,V[K-1]: 前K大奇异值对应的右奇异向量
*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include <stdbool.h>
#include "svd.h"

const int MAX_ITER = 100000;
const double eps = 0.0000001;

double get_norm(double *x, int n) {
	double r = 0;
	for (int i = 0;i<n;i++)
		r += x[i] * x[i];
	return sqrt(r);
}
double normalize(double *x, int n) {
	double r = get_norm(x, n);
	if (r<eps)
		return 0;
	for (int i = 0;i<n;i++)
		x[i] /= r;
	return r;
}

double product(double*a, double *b, int n) {
	double r = 0;
	for (int i = 0;i<n;i++)
		r += a[i] * b[i];
	return r;
}

void orth(double *a, double *b, int n) {
	double r = product(a, b, n);
	for (int i = 0;i<n;i++)
		b[i] -= r*a[i];

}

bool svd(double A[][N], int K, double U[][N], double S[N], double V[][N]) {
	srand(time(0));
	double *left_vector = (double*)malloc(sizeof(double) * M);
	double *next_left_vector = (double*)malloc(sizeof(double) * M);
	double *right_vector = (double*)malloc(sizeof(double) * N);
	double *next_right_vector = (double*)malloc(sizeof(double) * N);
	int col = 0;
	for (int col = 0;col<K;col++) {
		double diff = 1;
		double r = -1;
		while (1) {
			for (int i = 0;i<M;i++)
				left_vector[i] = (float)rand() / RAND_MAX;
			if (normalize(left_vector, M)>eps)
				break;
		}

		for (int iter = 0;diff >= eps && iter<MAX_ITER;iter++) {
			memset(next_left_vector, 0, sizeof(double)*M);
			memset(next_right_vector, 0, sizeof(double)*N);
			for (int i = 0;i<M;i++)
				for (int j = 0;j<N;j++)
					next_right_vector[j] += left_vector[i] * A[i][j];

			r = normalize(next_right_vector, N);
			if (r<eps) break;
			for (int i = 0;i<col;i++)
				orth(&V[i][0], next_right_vector, N);
			normalize(next_right_vector, N);

			for (int i = 0;i<M;i++)
				for (int j = 0;j<N;j++)
					next_left_vector[i] += next_right_vector[j] * A[i][j];
			r = normalize(next_left_vector, M);
			if (r<eps) break;
			for (int i = 0;i<col;i++)
				orth(&U[i][0], next_left_vector, M);
			normalize(next_left_vector, M);
			diff = 0;
			for (int i = 0;i<M;i++) {
				double d = next_left_vector[i] - left_vector[i];
				diff += d*d;
			}

			memcpy(left_vector, next_left_vector, sizeof(double)*M);
			memcpy(right_vector, next_right_vector, sizeof(double)*N);
		}
		if (r >= eps) {
			S[col] = r;
			memcpy((char *)&U[col][0], left_vector, sizeof(double)*M);
			memcpy((char *)&V[col][0], right_vector, sizeof(double)*N);
		}
		else {
			//cout << r << endl;
			break;
		}
	}
	free(next_left_vector);
	free(next_right_vector);
	free(left_vector);
	free(right_vector);

	return true;
}


/* Numerical Recipes standard error handler */
void nrerror(char error_text[])
{
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	//system("pause");
	exit(1);
}

/* allocate a double vector with subscript range v[nl..nh] */
double *fvector(long nl, long nh)
{
	double *v;
	v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof (double)));
	if (!v) nrerror("allocation failure in vector()");
	return v - nl + NR_END;
}

/* free a double vector allocated with vector() */
void free_fvector(double *v, long nl, long nh)
{
	free((FREE_ARG)(v + nl - NR_END));
}

// using sqrt instead of sqrt
double pythag(double a, double b)
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa*sqrt(1.0 + SQR(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + SQR(absa / absb)));
}

/* svd operation api */
void svdcmp(double **a, int m, int n, double w[], double **v)
{
	double pythag(double a, double b);
	int flag, i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z, *rv1;
	rv1 = fvector(1, n);

	g = scale = anorm = 0.0;
	for (i = 1; i <= n; i++)
	{
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m)
		{
			for (k = i; k <= m; k++) scale += fabs(a[k][i]);
			if (scale)
			{
				for (k = i; k <= m; k++)
				{
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++)
				{
					for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
				}
				for (k = i; k <= m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n)
		{
			for (k = l; k <= n; k++) scale += fabs(a[i][k]);
			if (scale)
			{
				for (k = l; k <= n; k++)
				{
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++)
				{
					for (s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++) a[j][k] += s * rv1[k];
				}
				for (k = l; k <= n; k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n; i >= 1; i--)
	{
		if (i < n)
		{
			if (g)
			{
				for (j = l; j <= n; j++)
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++)
				{
					for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++) v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = IMIN(m, n); i >= 1; i--)
	{
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++) a[i][j] = 0.0;
		if (g)
		{
			g = 1.0 / g;
			for (j = l; j <= n; j++)
			{
				for (s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
				f = (s / a[i][i]) * g;
				for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
			}
			for (j = i; j <= m; j++) a[j][i] *= g;
		}
		else for (j = i; j <= m; j++) a[j][i] = 0.0;
		++a[i][i];
	}
	for (k = n; k >= 1; k--)
	{
		for (its = 1; its <= 30; its++)
		{
			flag = 1;
			for (l = k; l >= 1; l--)
			{
				nm = l - 1;
				if ((double)(fabs(rv1[l]) + anorm) == anorm)
				{
					flag = 0;
					break;
				}
				if ((double)(fabs(w[nm]) + anorm) == anorm) break;
			}
			if (flag)
			{
				c = 0.0;

				s = 1.0;
				for (i = l; i <= k; i++)
				{
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((double)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++)
					{
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}
			z = w[k];
			if (l == k)
			{
				if (z < 0.0)
				{
					w[k] = -z;
					for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 100) nrerror("no convergence in 30 svdcmp iterations");
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++)
				{
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++)
				{
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_fvector(rv1, 1, n);
}


void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
	int jj, j, i;
	double s, *tmp;
	tmp = fvector(1, n);
	for (j = 1; j <= n; j++)
	{
		s = 0.0;
		if (w[j])
		{
			for (i = 1; i <= m; i++) s += u[i][j] * b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (j = 1; j <= n; j++)
	{
		s = 0.0;
		for (jj = 1; jj <= n; jj++) s += v[j][jj] * tmp[jj];
		x[j] = s;
	}
	free_fvector(tmp, 1, n);
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **fmatrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;
	/* allocate pointers to rows */
	m = (double **)malloc((size_t)((nrow + NR_END) * sizeof (double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof (double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

/* free a double matrix allocated by matrix() */
void free_fmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

