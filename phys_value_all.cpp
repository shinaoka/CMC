#include <stdlib.h>
#include "mpi.h"
#include "mfmemory.h"
#include "phys_value_all.h"

void init_phys_value_all(Phys_value_all *pva, int num_phys_value)
{
	int mfint[7]; //for memory allocator
	d_malloc2(pva->phys_value, num_phys_value, 2)
	d_malloc2(pva->phys_value_suscpt, num_phys_value, 2)
	d_malloc1(pva->E, 2)
	d_malloc1(pva->C, 2)
	d_malloc2(pva->SG_op, 2, 2)
}

void free_phys_value_all(Phys_value_all *pva, int num_phys_value)
{
	int mfint[7]; //for memory allocator
	d_free2(pva->phys_value, num_phys_value, 2)
	d_free2(pva->phys_value_suscpt, num_phys_value, 2)
	d_free1(pva->E, 2)
	d_free1(pva->C, 2)
	d_free2(pva->SG_op, 2, 2)
}

void cp_phys_value_all(Phys_value_all *src, Phys_value_all *dst, int num_phys_value)
{
	int i,j;

	for(i=0; i<=num_phys_value-1; i++)
	{
		for(j=0; j<=1; j++)
		{
			dst->phys_value[i][j] = src->phys_value[i][j];
			dst->phys_value_suscpt[i][j] = src->phys_value_suscpt[i][j];
		}
	}
	for(i=0; i<=1; i++)
	{
		dst->SG_op[i][0] = src->SG_op[i][0];
		dst->SG_op[i][1] = src->SG_op[i][1];
	}
	for(j=0; j<=1; j++)
	{
		dst->E[j] = src->E[j];
		dst->C[j] = src->C[j];
	}
}

#ifdef _mpi_use
/****
int MPI_Gather_phys_value_all(int i_beta, Phys_value_all *pva, Phys_value_all* t_pva, int num_phys_value, int num_thread, int my_rank)
{
	double *buffer, *rsv_buffer, *p_buffer;
	int len_buffer;
	int i, j, i_thread;
	Phys_value_all *t_pva2;

	len_buffer = 4*num_phys_value+5+4;
	d_malloc1(buffer, len_buffer)
	d_malloc1(rsv_buffer, num_thread*len_buffer)

	p_buffer = buffer;
	*p_buffer = (double) i_beta;
	p_buffer++;
	for(i=0; i<=num_phys_value-1; i++)
	{
		*p_buffer = t_pva->phys_value[i][0];
		p_buffer++;
		*p_buffer = t_pva->phys_value[i][1];
		p_buffer++;
		*p_buffer = t_pva->phys_value_suscpt[i][0];
		p_buffer++;
		*p_buffer = t_pva->phys_value_suscpt[i][1];
		p_buffer++;
	}
	*p_buffer = t_pva->E[0];
	p_buffer++;
	*p_buffer = t_pva->E[1];
	p_buffer++;
	*p_buffer = t_pva->C[0];
	p_buffer++;
	*p_buffer = t_pva->C[1];
	p_buffer++;

	for(i=0; i<=1; i++)
	{
		for(j=0; j<=1; j++)
		{
			*p_buffer = t_pva->SG_op[i][j];
			p_buffer++;
		}
	}

	MPI_Gather(buffer, len_buffer, MPI_DOUBLE, rsv_buffer, len_buffer, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(my_rank==0)
	{
		p_buffer = rsv_buffer;
		for(i_thread=0; i_thread<=num_thread-1; i_thread++)
		{
			t_pva2 = &pva[(int) *p_buffer];
			p_buffer++;
			for(i=0; i<=num_phys_value-1; i++)
			{
				t_pva2->phys_value[i][0] = *p_buffer;
				p_buffer++;
				t_pva2->phys_value[i][1] = *p_buffer;
				p_buffer++;
				t_pva2->phys_value_suscpt[i][0] = *p_buffer;
				p_buffer++;
				t_pva2->phys_value_suscpt[i][1] = *p_buffer;
				p_buffer++;
			};
			t_pva2->E[0] = *p_buffer;
			p_buffer++;
			t_pva2->E[1] = *p_buffer;
			p_buffer++;
			t_pva2->C[0] = *p_buffer;
			p_buffer++;
			t_pva2->C[1] = *p_buffer;
			p_buffer++;
			for(i=0; i<=1; i++)
			{
				for(j=0; j<=1; j++)
				{
					t_pva2->SG_op[i][j] = *p_buffer;
					p_buffer++;
				}
			}
		}
	}

	d_free1(buffer, len_buffer)
	d_free1(rsv_buffer, num_buffer*len_buffer)
	return 0;
}
****/
#endif
