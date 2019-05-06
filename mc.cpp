#include "mc.h"
#ifdef _mpi_use
	#include "mpi.h"
	#include "mpi_util.h"
#endif

#include <iostream>
#include <fstream>
#include <sstream>

#undef DEBUG

#include "exchange.cpp"
#include "spin.h"

using namespace std;

/*-------------------------------------------------------------
 * MC code for Heisenberg spin models
 * phys_value: [0][0:1]    <Mx>, error of <Mx>
 *             [1][0:1]    <My>, error of <My>
 *             [2][0:1]    <Mz>, error of <Mz>
 *             [3][0:1]    <Q_{3z^2-r-2}>
 *             [4][0:1]    <Q_{x^2-y2}>
 *             [5][0:1]    <Q_{xy}>
 *             [6][0:1]    <Q_{xz}>
 *             [7][0:1]    <Q_{yz}>
 *             [8][0:1]    <Q^2>
 *
 * time (sec): [0] overlap cluster
 *             [1] signle-spin flip
 *             [2] loop move
 *             [3] phys value
 *             [4] exchange
 * phys_value_suscpt: susceptibility and error of each correspoinding physical quantity
 *-------------------------------------------------------------*/
int mc(int i_def, int i_run, dsfmt_t *dsfmt, Ham *ham, int num_beta, double *beta_list,
		int num_sbl_mag, vector<complex<double> > sbl_mag_coeff, 
                std::vector<CROSS_PRODUCT_TYPE>& cross_product_coeff, 
		int num_replica, int num_seperated_replica_set, int num_MC_step,
		int num_MC_step_thermalization, int MC_step_per_ave, int MC_step_per_ave_ss,
		MCUpdateParam *mc_update_param, 
		double **ex_rate, double *q2_dynamics, vector<PhysValueRecorder>& pvr, vector<PhysValueRecorder>& sr,
		int num_ss_pair, int **ss_pair, vector<PhysValueRecorder>& sscr, 
 		vector<vector<double> >& auto_correlation, vector<double> &time, bool abs_phys)
{
	const int num_pow_phys = 4;

	/*variable declaration-------------------------------------*/
	double beta;
	int i_site, i_site2, i_spin, i_coordinate, i_beta; //loop variables
	int i, j, k, itmp, itmp2, itmp3, itmp4, itmp5; //general-purpose temporary variable
	int Ns; // Number of sites
	double rtmp1, rtmp2, rtmp3; //general-purpose temporary variable
	int mfint[7]; //for memory allocator
	double *config;
	int num_beta_per_process;
	GraphNetwork graph_n;

	Replica *replica, *p_replica;
	int num_independent_replica_pair, num_replica_per_set, i_replica,
			i_replica2, index_replica_set;

	double t_phys_value[NUM_PHYS_VALUE];
	PhysValueCalc **phys_value_calc_list, **SG_calc_list, **SSC_calc_list;

	/*Ising anisotropy-----------------------------------------*/
	//double **ia;

#ifdef _utime
	//for usage time
	double ttime[100];
#endif

	//MPI
	int my_rank, p_num_check, len_buffer_ex;
	double *send_buffer_ex, *recv_buffer_ex;

	//Replica exchange
	int *ex_counter;
	int INTV_BETA_UPDATE;

	// update statics
	//0: single-spin flip
	//1: over relaxation
	//2: spin inversion
	//3: loop_move num_site_visited
	//4: loop_move num_loop_found
	//5: loop_move num_loop_flipped
	//6: loop_move num_site_flipped
	//7: overlap cluster
	unsigned long **update_counter;

	/**** work arrays for loop move              ****/
	LoopMoveWorkArrays lm_wrk;

	/**** work arrays for auto correlations ****/
	double ***ref_config;

	//for MC loop
	double d_config[3], phi, cos_theta, dE;

	const double PI = 3.1415926535;

	/*calc parameters------------------------------------------*/
#ifdef _mpi_use
	MPI_Comm_rank(MY_MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MY_MPI_COMM_WORLD, &p_num_check);
#else
	my_rank = 0;
	p_num_check = 1;
#endif
	Ns = ham->Ns;
	num_beta_per_process = num_beta / p_num_check;
	INTV_BETA_UPDATE = 400;
	num_replica_per_set = num_replica / num_seperated_replica_set;
	num_independent_replica_pair = num_replica_per_set
			* num_seperated_replica_set * (num_seperated_replica_set - 1) / 2;
	int num_sampling = num_MC_step/MC_step_per_ave;
	int num_sampling_ss = num_MC_step/MC_step_per_ave_ss;

	/*check consistency of input parameters--------------------*/
	if (num_MC_step%MC_step_per_ave != 0 || num_MC_step_thermalization%MC_step_per_ave != 0)
		return -1;
	if (num_MC_step%MC_step_per_ave_ss != 0 || num_MC_step_thermalization%MC_step_per_ave_ss != 0)
		return -1;
	
	/*init phys_value_recorder---------------------------------*/
	replica = (Replica*) malloc(num_beta_per_process * num_replica
			* sizeof(Replica));
	{
		phys_value_calc_list = new PhysValueCalc*[num_beta_per_process];
		SG_calc_list = new PhysValueCalc*[num_beta_per_process];
		SSC_calc_list = new PhysValueCalc*[num_beta_per_process];

		for(int i_beta = 0; i_beta < num_beta_per_process; ++ i_beta)
		{
			phys_value_calc_list[i_beta] = new BuiltInPhysValueCalc(num_replica, ham,
                                num_sbl_mag, sbl_mag_coeff,
                                cross_product_coeff,
                                abs_phys);
			SG_calc_list[i_beta] = new SGCalc(num_replica, num_seperated_replica_set, ham);
			SSC_calc_list[i_beta] = new SSCCalc(num_replica, ham, num_ss_pair, ss_pair);

			for(int i_replica = 0; i_replica < num_replica; ++ i_replica)
			{
				phys_value_calc_list[i_beta]->add_ptr_replica(i_replica, &replica[i_beta + i_replica * num_beta_per_process]);
				SG_calc_list[i_beta]->add_ptr_replica(i_replica, &replica[i_beta + i_replica * num_beta_per_process]);
				SSC_calc_list[i_beta]->add_ptr_replica(i_replica, &replica[i_beta + i_replica * num_beta_per_process]);
			}
			int num_sample_disregarded = num_MC_step_thermalization/MC_step_per_ave;
			PhysValueRecorder *ptr = new PhysValueRecorder(phys_value_calc_list[i_beta], num_pow_phys, num_sampling, num_sample_disregarded, 1);
			pvr[i_beta] = *ptr;
			delete ptr;

			ptr = new PhysValueRecorder(SG_calc_list[i_beta], 2, num_sampling, num_sample_disregarded, 1);
			sr[i_beta] = *ptr;
			delete ptr;

			int num_sample_disregarded_ss = num_MC_step_thermalization/MC_step_per_ave_ss;
			ptr = new PhysValueRecorder(SSC_calc_list[i_beta], 3, num_sampling_ss, num_sample_disregarded_ss, 0);
			sscr[i_beta] = *ptr;
			delete ptr;
		}
	}

	/*create work spaces for auto correlations-----------------*/
	auto_correlation.resize(num_beta_per_process);
	for(int i_beta = 0; i_beta < num_beta_per_process; ++ i_beta)
		auto_correlation[i_beta].resize(num_MC_step-num_MC_step_thermalization);
	malloc3_seq(ref_config, num_replica*num_beta_per_process, ham->Ns, 3);

	/****
	d_malloc2(ia, Ns, 4)
	//default easy axis is z axis
	for (i_site = 0; i_site <= Ns - 1; i_site++){
		ia[i_site][0] = 0.0;

		ia[i_site][1] = 0.0;
		ia[i_site][2] = 0.0;
		ia[i_site][3] = 1.0;
	}
	//overwrite default easy axis
	for (i = 0; i <= ham->num_ia - 1; i++){
		ia[ham->ia_site[i]][0] = ham->ia[i][0]; //alpha
		ia[ham->ia_site[i]][1] = ham->ia[i][1]; //r_x
		ia[ham->ia_site[i]][2] = ham->ia[i][2]; //r_y
		ia[ham->ia_site[i]][3] = ham->ia[i][3]; //r_z
	}
	//normalize
	for (i_site = 0; i_site <= Ns - 1; i_site++){
		double rtmp = 1.0/sqrt(norm2(&ia[i_site][1]));
		for(int mu = 0; mu < 3; ++ mu)
			ia[i_site][mu+1] *= rtmp;
	}
	****/

	/*devide lattice into a set of complete graph----------------------------*/
	if(mc_update_param->max_num_site_loop > 0)
	{
		construct_complete_graph_network(&graph_n, dsfmt, ham,
			ham->idx_cnct_bond, ham->idx_cnctd_site, ham->coord_num);
		printf("%d complete graph found\n", graph_n.num_graph);
		allocate_loop_move_work_arrays(&lm_wrk, Ns, graph_n.num_graph, graph_n.max_num_nn_graph);
	}

	/*initialize work array--------------------------------------------------*/
#ifdef DEBUG
	printf("my_rank=%d initializing work arrays...\n", my_rank);
#endif
#ifdef _utime
	time[0] = 0.0;
	time[1] = 0.0;
	time[2] = 0.0;
	time[3] = 0.0;
	time[4] = 0.0;
#endif
	i_malloc1(ex_counter, num_beta_per_process)
	malloc2_seq(update_counter, num_beta_per_process, 9); //0, 1, 2=>single spin flip, 3=>loop move, 4=>over-lap cluster
	if (p_num_check > 1)
	{
		len_buffer_ex = 3 * Ns + 2;
		d_malloc1(send_buffer_ex, len_buffer_ex);
		d_malloc1(recv_buffer_ex, len_buffer_ex);
	}
	for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
	{
		ex_counter[i_beta] = 0;
		set_val1(update_counter[i_beta], 0UL, 9);
	}

	/*initial configuration-------------------------------------*/
	for (i_replica = 0; i_replica <= num_replica - 1; i_replica++)
	{
		for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
		{
			p_replica = &replica[i_beta + i_replica * num_beta_per_process];
			d_malloc2(p_replica->config, Ns, 3);
		}
	}

	if(mc_update_param->load_config == 1){
		ostringstream oss;
		oss << "config-def" << i_def << "-run" << i_run << "-p" << my_rank << ".txt";
		FILE *fp = fopen(oss.str().c_str(), "r");
		cout << "Loading config..." << endl;
		for (i_replica = 0; i_replica <= num_replica - 1; i_replica++)
		{
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
			{
				p_replica = &replica[i_beta + i_replica * num_beta_per_process];
				double **p_config = replica[i_beta + i_replica * num_beta_per_process].config;
				for (i_site = 0; i_site < Ns; ++ i_site)
				{
					fscanf(fp, "%lf %lf %lf",
						&p_config[i_site][0], &p_config[i_site][1], &p_config[i_site][2]);
				}
				p_replica->E = calc_E(ham, p_config);
				p_replica->id = i_beta + my_rank * num_beta_per_process;
			}
		}
		fclose(fp);
	} else if (mc_update_param->load_config == 1) {
		for (i_replica = 0; i_replica <= num_replica - 1; i_replica++){
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++){
				p_replica = &replica[i_beta + i_replica * num_beta_per_process];
				for (i_site = 0; i_site < Ns; ++ i_site){
					if (ham->spin[i_site] == ISING_SPIN){
						double inv_norm2 = 1.0/sqrt(norm2(ham->ia[i_site]+1));
						double dsign = 2.0 * (dsfmt_genrand_uint32(dsfmt) % 2) - 1.0;
						for(int mu = 0; mu < 3; ++mu)
							p_replica->config[i_site][mu] = inv_norm2*dsign*ham->ia[i_site][mu+1];
					} else if (ham->spin[i_site] == HEISENBERG_SPIN) {
						phi = 2.0 * PI * dsfmt_genrand_close_open(dsfmt);
						cos_theta = dsfmt_genrand_close_open(dsfmt) * 2.0 - 1.0;
						p_replica->config[i_site][0] = sqrt(1.0 - pow(cos_theta,
								2.0)) * cos(phi);
						p_replica->config[i_site][1] = sqrt(1.0 - pow(cos_theta,
								2.0)) * sin(phi);
						p_replica->config[i_site][2] = cos_theta;
					} else if (ham->spin[i_site] == XY_SPIN) {
						double theta = 2.0 * PI * dsfmt_genrand_close_open(dsfmt);
						p_replica->config[i_site][0] = std::cos(theta);
						p_replica->config[i_site][1] = std::sin(theta);
						p_replica->config[i_site][2] = 0.0;
					}
				}
				p_replica->E = calc_E(ham, p_replica->config);
				p_replica->id = i_beta + my_rank * num_beta_per_process;
			}
		}
	} else {
		throw std::runtime_error("Invalid value of mc_update_param->load_config");
	}

	/*MC loop--------------------------------------------------*/
#ifdef DEBUG
	printf("my_rank=%d starting MC loops...\n", my_rank);
#endif
	for (int i_MC_step = 0; i_MC_step <= num_MC_step - 1; ++ i_MC_step)
	{
#ifdef DEBUG
		printf("my_rank=%d i_MC_step=%d\n", my_rank, i_MC_step);
#endif

#ifdef _utime
		ttime[0] = getrusage_sec();
#endif
		if(i_MC_step == num_MC_step_thermalization)
		{
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; ++ i_beta)
				ex_counter[i_beta] = 0;
		}

		/**** optimize beta ****/
		if(mc_update_param->num_ex_opt > 0)
		{
			if(i_MC_step > 0 && i_MC_step%mc_update_param->ex_opt_interval == 0 && i_MC_step/mc_update_param->ex_opt_interval <= mc_update_param->num_ex_opt)
			{
				opt_beta(num_beta, num_beta_per_process, my_rank, beta_list
					, mc_update_param->ex_opt_min_beta_index, mc_update_param->ex_opt_max_beta_index, ex_counter);
				for (i_beta = 0; i_beta <= num_beta_per_process - 1; ++ i_beta)
					ex_counter[i_beta] = 0;
			}
		}
		
#ifdef _utime
		ttime[1] = getrusage_sec();
		time[0] += ttime[1] - ttime[0];
#endif

		/**** START OF LOCAL UPDATE ****/
		for (i_replica = 0; i_replica <= num_replica - 1; i_replica++)
		{
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
			{
#ifdef _utime
				ttime[0] = getrusage_sec();
#endif
				beta = beta_list[i_beta + my_rank * num_beta_per_process];
				p_replica = &replica[i_beta + i_replica * num_beta_per_process];

				if(beta <= mc_update_param->max_beta_single_spin_flip)
				{
					update(dsfmt, p_replica, beta, ham, ham->idx_cnct_bond,
							ham->idx_cnctd_site, ham->coord_num,
							ham->ia, ham->mag, mc_update_param->num_overrelaxation,
							update_counter[i_beta], mc_update_param->fluctuation_random_walk, false);
					/****
					if(fabs(calc_E(ham, p_replica->config)-p_replica->E)>1E-8){
						cout << calc_E(ham, p_replica->config) << p_replica->E << endl;
					}
					****/
				}

#ifdef _utime
				ttime[1] = getrusage_sec();
#endif

				if(beta >= mc_update_param->min_beta_loop_move && mc_update_param->max_num_site_loop > 0)
				{
					if(mc_update_param->max_num_site_loop > 0 && graph_n.num_graph>0){
						loop_move_ISING(dsfmt, p_replica, beta, ham, ham->idx_cnct_bond,
							ham->idx_cnctd_site, ham->coord_num, ham->ia,
							mc_update_param->max_num_site_loop, mc_update_param->num_update_fix_loop,
							mc_update_param->num_site_ising_axis, mc_update_param->loop_flip_method,
							&graph_n, &lm_wrk, &itmp, &itmp2, &itmp3, &itmp4, &itmp5);
						update_counter[i_beta][3] += itmp;
						update_counter[i_beta][4] += itmp2;
						update_counter[i_beta][5] += itmp3;
						update_counter[i_beta][6] += itmp4;
						update_counter[i_beta][7] += itmp5;
					}
				}

#ifdef _utime
				ttime[2] = getrusage_sec();
				time[1] += ttime[1] - ttime[0];
				time[2] += ttime[2] - ttime[1];
#endif
			}
		}
		/**** END OF LOCAL UPDATE ****/

		/**** sampling physical quantities ****/
#ifdef _utime
		ttime[0] = getrusage_sec();
#endif
		/**** auto correlation ****/
		if(mc_update_param->auto_correlation != 0)
		{
			if(i_MC_step - num_MC_step_thermalization == 0)
			{
				for(int i_replica = 0; i_replica < num_replica*num_beta_per_process; ++ i_replica)
				{
					for(int i_site = 0; i_site < ham->Ns; ++ i_site)
					{
						for(int ix = 0; ix < 3; ++ ix)
							ref_config[i_replica][i_site][ix] = replica[i_replica].config[i_site][ix];
					}
				}
			}
			if(i_MC_step - num_MC_step_thermalization >= 0)
			{
				for(int i_beta = 0; i_beta < num_beta_per_process; ++ i_beta)
				{
					auto_correlation[i_beta][i_MC_step - num_MC_step_thermalization] = 0.0;
					for(int i_replica = 0; i_replica < num_replica; ++ i_replica)
					{
						double t_overlap = 0.0;
						for(int i_site = 0; i_site < ham->Ns; ++ i_site)
						{
							for(int ix = 0; ix < 3; ++ ix)
							{
								t_overlap += ref_config[i_beta + i_replica*num_beta_per_process][i_site][ix]
									*replica[i_beta + i_replica*num_beta_per_process].config[i_site][ix];
							}
						}
						auto_correlation[i_beta][i_MC_step - num_MC_step_thermalization] += t_overlap*t_overlap;
					}
					auto_correlation[i_beta][i_MC_step - num_MC_step_thermalization] /= (1.0*ham->Ns*num_replica);
				}
			}
		}

		if(i_MC_step%MC_step_per_ave == 0){	
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++){
				pvr[i_beta].sampling();
				sr[i_beta].sampling();
			}
		}
		if(i_MC_step%MC_step_per_ave_ss == 0){	
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
				sscr[i_beta].sampling();
		}
#ifdef _utime
		ttime[1] = getrusage_sec();
		time[3] += ttime[1] - ttime[0];
#endif

		/**** START OF EXCHANGE ****/
#ifdef _utime
		ttime[0] = getrusage_sec();
#endif
		if (mc_update_param -> ex_interval > 0)
		{
			if (num_beta > 1 && i_MC_step % mc_update_param->ex_interval == 0)
			{
				for (i_replica = 0; i_replica <= num_replica - 1; i_replica++)
				{
					exchange(dsfmt, Ns, num_beta_per_process,
							&replica[i_replica * num_beta_per_process],
							p_num_check, my_rank, num_beta, beta_list,
							mc_update_param->num_inner_process_ex, ex_counter, send_buffer_ex, recv_buffer_ex);
				}
			}
		}
#ifdef _utime
		ttime[1] = getrusage_sec();
		time[4] += ttime[1] - ttime[0];
#endif
		/**** END OF EXCHANGE ****/
	}
	/*end MC loop----------------------------------------------*/

#ifdef DEBUG
	printf("my_rank=%d avearaging physical quantities...\n", my_rank);
#endif

	for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
	{
		beta = beta_list[i_beta + my_rank * num_beta_per_process];
		/*exchange exchange rate----------------------------*/
		ex_rate[i_beta][0] = ((double) ex_counter[i_beta]) / (num_replica*(num_MC_step-num_MC_step_thermalization)/mc_update_param->ex_interval);
		ex_rate[i_beta][1] = ((double) update_counter[i_beta][0])
				/ (1.0 * num_replica * num_MC_step * Ns);
		ex_rate[i_beta][2] = ((double) update_counter[i_beta][1])
				/ (1.0 * num_replica * num_MC_step * Ns);
		ex_rate[i_beta][3] = ((double) update_counter[i_beta][2])
				/ (1.0 * num_replica * num_MC_step * Ns);
		ex_rate[i_beta][4] = ((double) update_counter[i_beta][3])
				/ (num_replica * num_MC_step);
		ex_rate[i_beta][5] = ((double) update_counter[i_beta][4])
				/ (num_replica * num_MC_step);
		ex_rate[i_beta][6] = ((double) update_counter[i_beta][5])
				/ (num_replica * num_MC_step);
		ex_rate[i_beta][7] = ((double) update_counter[i_beta][6])
				/ (num_replica * num_MC_step);
		ex_rate[i_beta][8] = ((double) update_counter[i_beta][7])
				/ (num_replica * num_MC_step);
		ex_rate[i_beta][9] = ((double) update_counter[i_beta][8])
				/ (num_seperated_replica_set * num_replica_per_set * num_MC_step);
	}

	/*dump spin configurations------------------------------------------------------------------------*/
	if(mc_update_param->dump_config == 1)
	{
		ostringstream oss;
		oss << "config-def" << i_def << "-run" << i_run << "-p" << my_rank << ".txt";
		FILE *fp = fopen(oss.str().c_str(), "w");
		cout << "Dumping config..." << endl;
		for (i_replica = 0; i_replica <= num_replica - 1; i_replica++)
		{
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
			{
				double **p_config = replica[i_beta + i_replica * num_beta_per_process].config;
				for (i_site = 0; i_site < Ns; ++ i_site)
				{
					fprintf(fp, "%25.16le %25.16le %25.16le\n",
						p_config[i_site][0], p_config[i_site][1], p_config[i_site][2]);
				}
			}
		}
		fclose(fp);
	}

	/*free work arrays------------------------------------------------------------------------*/
	i_free1(ex_counter, num_beta_per_process)
	for (i_replica = 0; i_replica <= num_replica * num_beta_per_process - 1; i_replica++)
	{
		d_free2(replica[i_replica].config, Ns, 3);
	}
	free(replica);
	free2_seq(update_counter);
	if (p_num_check > 1)
	{
		d_free1(send_buffer_ex, len_buffer_ex);
		d_free1(recv_buffer_ex, len_buffer_ex);
	}

	/*free work spaces for auto correlations-----------------*/
	free3_seq(ref_config);

	if(mc_update_param->max_num_site_loop > 0 && graph_n.num_graph>0)
	{
		free_complete_graph_network(&graph_n);
		free_loop_move_work_arrays(&lm_wrk, Ns, graph_n.num_graph, graph_n.max_num_nn_graph);
	}

	delete[] phys_value_calc_list;

	return 0;
}
