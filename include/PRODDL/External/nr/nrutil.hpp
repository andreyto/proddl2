#ifndef AT_NR_UTILS_H_

#define AT_NR_UTILS_H_

#include <stdlib.h>

namespace PRODDL { namespace nr {


  template<typename T_num> T_num SQR(T_num a) { a*a; }
  template<typename T_num> T_num DSQR(T_num a) { a*a; }

template<typename T_num> T_num DMAX(T_num a, T_num b) { return a > b ? a : b; }

template<typename T_num> T_num DMIN(T_num a, T_num b) { return a < b ? a : b; }

template<typename T_num> T_num FMAX(T_num a, T_num b) { return a > b ? a : b; }

template<typename T_num> T_num FMIN(T_num a, T_num b) { return a < b ? a : b; }

template<typename T_num> T_num LMAX(T_num a, T_num b) { return a > b ? a : b; }

template<typename T_num> T_num LMIN(T_num a, T_num b) { return a < b ? a : b; }

template<typename T_num> T_num IMAX(T_num a, T_num b) { return a > b ? a : b; }

template<typename T_num> T_num IMIN(T_num a, T_num b) { return a < b ? a : b; }

template<typename T_num> T_num SIGN(T_num a, T_num b) { return b >= 0 ? fabs(a) : -fabs(b); }


 typedef char* FREE_ARG;
  enum { NR_END = 1 };

void nrerror(char error_text[]);

int *ivector(long nl, long nh);

unsigned char *cvector(long nl, long nh);

unsigned long *lvector(long nl, long nh);

double *dvector(long nl, long nh);

double **dmatrix(long nrl, long nrh, long ncl, long nch);

int **imatrix(long nrl, long nrh, long ncl, long nch);

void free_ivector(int *v, long nl, long nh);

void free_cvector(unsigned char *v, long nl, long nh);

void free_lvector(unsigned long *v, long nl, long nh);

void free_dvector(double *v, long nl, long nh);

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);


// template functions

template<typename T_num>
T_num *vector(long nl, long nh)
/* allocate a T_num vector with subscript range v[nl..nh] */
{
	T_num *v;

	v=(T_num *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(T_num)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

template<typename T_num>
T_num **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a T_num matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	T_num **m;

	/* allocate pointers to rows */
	m=(T_num **) malloc((size_t)((nrow+NR_END)*sizeof(T_num*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(T_num *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(T_num)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

template<typename T_num>
T_num **submatrix(T_num **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	T_num **m;

	/* allocate array of pointers to rows */
	m=(T_num **) malloc((size_t) ((nrow+NR_END)*sizeof(T_num*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

template<typename T_num>
T_num **convert_matrix(T_num *a, long nrl, long nrh, long ncl, long nch)
/* allocate a T_num matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	T_num **m;

	/* allocate pointers to rows */
	m=(T_num **) malloc((size_t) ((nrow+NR_END)*sizeof(T_num*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

template<typename T_num>
T_num ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a T_num 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	T_num ***t;

	/* allocate pointers to pointers to rows */
	t=(T_num ***) malloc((size_t)((nrow+NR_END)*sizeof(T_num**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(T_num **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(T_num*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(T_num *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(T_num)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

template<typename T_num>
void free_vector(T_num *v, long nl, long nh)
/* free a T_num vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

template<typename T_num>
void free_matrix(T_num **m, long nrl, long nrh, long ncl, long nch)
/* free a T_num matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

template<typename T_num>
void free_submatrix(T_num **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

template<typename T_num>
void free_convert_matrix(T_num **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

template<typename T_num>
void free_f3tensor(T_num ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a T_num f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

} } // namespace PRODDL::nr

#endif /* AT_NR_UTILS_H_ */

