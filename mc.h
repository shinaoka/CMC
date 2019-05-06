/*-------------------------------------------------------------
 * MC code for Heisenberg spin models
 * phys_value: [0][0,1] <Mx>, error of <Mx>
 *             [1][0,1] <My>, error of <My>
 *             [2][0,1] <Mz>, error of <Mz>
 *
 * phys_value_suscpt: susceptibilities and theirs errors for correspoinding physical quantities
 *-------------------------------------------------------------*/
#ifndef MC_H
#define MC_H

#ifdef _mpi_use
#include <mpi.h>
#endif

#include <math.h>
#include <vector>
#include <complex>
using namespace std;
#include "dSFMT.h"
#include "ham.h"
#include "replica.h"
#include "loop_move.h"
#include "phys_value_recorder.h"
#include "update.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

struct MCUpdateParam
{
	double max_beta_single_spin_flip, fluctuation_random_walk, fluctuation_OR;
	int num_overrelaxation;

  //for chaining
	int load_config, dump_config;

  //for parallel tempering
	int ex_interval, num_inner_process_ex;
  //for optimzation of parallel tempering
	int ex_opt_interval, num_ex_opt, ex_opt_min_beta_index, ex_opt_max_beta_index;
	//for loop move
	double min_beta_loop_move;
	int max_num_site_loop, num_update_fix_loop;
	//0: xyz
	//1: rotate
	//2: ising axis
	int loop_flip_method;
	int num_site_ising_axis;
	//Overlap-cluster update
	/*
	* 0: Swendsen-Wang algorithm for ISING (always accepted)
	* 1: Wolff algorithm (for ISING and Heisenberg) 
	*/
	int cluster_update_algorithm;
	int num_cluster_update;
	int auto_correlation;
};
typedef struct MCUpdateParam MCUpdateParam;


int mc(int i_def, int i_run, dsfmt_t *dsfmt, Ham *ham, int num_beta, double *beta_list, 
	int num_sbl_mag, vector<complex<double > > sbl_mag_coeff,
        std::vector<CROSS_PRODUCT_TYPE>& cross_product_coeff, 
        int num_replica, int num_seperated_replica_set,
	int num_MC_step, int num_MC_step_thermalization, int MC_step_per_ave, int MC_step_per_ave_ss, 
	MCUpdateParam *mc_update_param, double **ex_rate, double *q2_dynamics, vector<PhysValueRecorder>& pvr, vector<PhysValueRecorder>& sr,
	int num_ss_pair, int **ss_pair, vector<PhysValueRecorder>& ssr,
	vector<vector<double> >& auto_correlation,
	vector<double> &time, bool abs_phys);

void init_replica(Replica *rp, int Ns);
#endif
