/*-------------------------------------------------------------
 * util.h
 *-------------------------------------------------------------*/
#ifndef UTIL_H
#define UTIL_H

#include <string.h>
#include <vector>
using namespace std;

int mymod(int i, int N);
double dot_product(double *S0, double *S1);
double dmulti(double *darray, int start, int end, double C);
double dsum(double *darray, int start, int end);
double dsum_m(double *darray, double *darray2, int start, int end);
double getrusage_sec();

/****
inline void set_val_d1(vector<double>& array, double val)
{
	vector<double>::iterator it;

	while(it != array.end())
	{
		*it = val;
		++ it;
	}
}
****/

template<class T>
T cross_prod(T& a, T& b)
{
	T c(3); 
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
	return c;
}

inline double SQUARE(double x)
{
	return x * x;
}

inline double norm2(double *vec)
{
	return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
}

inline void set_val_d2(int **array, int val, int num_element1, int num_element2)
{
	int i, j;

	for(i=0; i <= num_element1 - 1 ;i ++)
	{
		for(j=0; j <= num_element2 - 1 ;j ++)
		{
			array[i][j] = val;
		}
	}
}

template<typename T>
inline void set_val1(T*& array, T val, int num_element)
{
	for(int i=0; i<= num_element -1 ;i ++)
	{
		*(array+i) = val;
	}
}


inline void set_val_i1(int *array, int val, int num_element)
{
	int i;

	for(i=0; i<= num_element -1 ;i ++)
	{
		*(array+i) = val;
	}
}

inline int isign(double rtmp)
{
	if(rtmp > 0.0)
	{
		return 1;
	}
	else if(rtmp < 0.0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}


#define icnv2(i, j, I, J) \
 (j)+(J)*(i)

#define icnv3(i, j, k, I, J, K) \
 (k)+(K)*(j)+(J)*(K)*(i)

#endif
