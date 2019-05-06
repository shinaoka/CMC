#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
/*-------------------------------------------------------------
 * util.c
 *-------------------------------------------------------------*/
int mymod(int i, int N)
{
	if(i>=0){
		return i%N;
	}else{
		return (i+N*((-i)/N+1))%N;
	}
}

double dot_product(double *S0, double *S1)
{
	return S0[0]*S1[0] + S0[1]*S1[1] + S0[2]*S1[2];
}

void dmulti(double *darray, int start, int end, double C)
{
	int i;

	for(i = start; i <= end; i++)
	{
		darray[i] = darray[i]*C;
	}
}

double dsum(double *darray, int start, int end)
{
	double tmp;
	int i;

	tmp = 0.0;
	for(i = start; i <= end; i++)
	{
		tmp = tmp+darray[i];
	};
	return tmp;
}

double dsum_m(double *darray, double *darray2, int start, int end)
{
	double tmp;
	int i;

	tmp = 0.0;
	for(i = start; i <= end; i++)
	{
		tmp = tmp+darray[i]*darray2[i];
	};
	return tmp;
}

double getrusage_sec()
{
	struct rusage t;
	struct timeval tv;
	getrusage(RUSAGE_SELF, &t);
	tv = t.ru_utime;
	return tv.tv_sec + (double)tv.tv_usec*1e-6;
}
