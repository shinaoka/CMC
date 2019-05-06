#ifndef PHYS_VALUE_ALL_H
#define PHYS_VALUE_ALL_H

#include <mpi.h>

typedef struct {
	double *E, *C, **phys_value, **phys_value_suscpt, **SG_op;
} Phys_value_all;

void init_phys_value_all(Phys_value_all *pva, int num_phys_value);
void free_phys_value_all(Phys_value_all *pva, int num_phys_value);
void cp_phys_value_all(Phys_value_all *src, Phys_value_all *dst, int num_phys_value);
#ifdef _mpi_use
int MPI_Gather_phys_value_all(int i_beta, Phys_value_all *pva, Phys_value_all* t_pva, int num_phys_value, int num_thread, int my_rank);
#endif
#endif
