#ifndef UPDATE_H
#define UPDATE_H

#include <stdlib.h>
#include <math.h>
#include "dSFMT.h"
#include "util.h"
#include "ham.h"
#include "replica.h"
#include "spin.h"

double calc_dE(double S0, double S1, double S2, double new_S0,
		double new_S1, double new_S2, double *local_h, double B00, double B01,
		double B02, double B11, double B12, double B22);

void calc_local_field(int i_site, Ham *ham, double **config,
		int *num_nn_two_body_int, int **idx_cnct_bond,
		int **idx_cnctd_site, double **ia, double *local_h, double *B00,
		double *B01, double *B02, double *B11, double *B12, double *B22);

void update(dsfmt_t *dsfmt, Replica *p_replica, double beta,
		Ham *ham, int **idx_cnct_bond,
		int **pointer_to_counterpart_Jij, int *coord_num,
		double **ia, double **mag, int num_OR, unsigned long *update_counter, double t_fluctuation_random_walk, bool dry_run);
#endif
