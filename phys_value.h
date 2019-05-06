#ifndef PHYS_VALUE_H
#define PHYS_VALUE_H

#define NUM_PHYS_VALUE 10
#include "util.h"
#include "replica.h"
#include "ham.h"
#include <math.h>
#include <vector>
#include <complex>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>


double calc_E(Ham *ham, double **config);
void calc_phys_value(int Ns, double **config, double *phys_value, int num_ss_pairs, int **ss_pairs,
	int num_sbl_mag, vector<complex<double > > &sbl_mag_coeff, vector<complex<double > > &sbl_mag, 
	vector<CROSS_PRODUCT_TYPE>& def_cross_product, vector<complex<double> >& cross_product, bool abs_phys);
void calc_squared_q(int Ns, Replica *replica1, Replica *replica2, double *q2, double *q4);
void average(int num_data, double *data, double *ave, double *sigma);

#endif
