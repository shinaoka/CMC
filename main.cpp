#ifdef _mpi_use
#include "mpi.h"
#include "mpi_util.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <complex>
#include "ham.h"
#include "mfmemory.h"
#include "mc.h"
#include "phys_value.h"
#include "mpi_util.h"
#include "phys_value_recorder.h"
#include <time.h>
#include <vector>
#include "dSFMT.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>


using namespace std;

#undef DEBUG

int main(int argc, char* argv[])
{
	Ham *pham;
	int Ns;
	double *beta, *beta_rw;
	int num_replica, num_seperated_replica_set, beta_mode, num_MC_step,
			num_MC_step_thermalization, MC_step_per_ave, MC_step_per_ave_ss, num_beta, num_beta_rw,
			num_beta_per_process;
	int calc_mode, dump_phys_value, num_run, index_start_run, output_phys_value_dynamics, num_pdef;
  int nine_comp_ex_mode;
	int i_site, i_coordinate, j, i_beta, itmp;
	unsigned long i_clock, i_seed;
	double **ex_rate;
	int mfint[7]; //for memory allocator
	FILE *iJp, *ipp;
	char t_char[1000];

	double rtmp, rtmp2, *ra_tmp, *ra_tmp2;

	/**** spin-spin correlations ****/
	int num_ss, **ss_pair;
	vector<PhysValueRecorder> ssr, t_ssr;

	vector<PhysValueRecorder> pvr, t_pvr;
	vector<PhysValueRecorder> sr, t_sr;

	vector<vector<double> > auto_correlation;
	double *t_q2_dynamics;

	/**** user-defined sublattice magnetization ****/
	int num_sbl_mag;
	vector<complex<double > > sbl_mag_coeff;


        /*** user-defined order parameters defined by cross product ***/
	vector<CROSS_PRODUCT_TYPE> cross_product_coeff;


	//Controll parameter for MC updates
	MCUpdateParam mc_update_param;

	//MPI
	int my_rank, p_num_check, my_world, my_splitted_rank, my_world_size;

	//Pseudo random generator
	dsfmt_t dsfmt;

	//usage time
	vector<double> utime;
	
	//abs for phys values
	int i_abs_phys;
	bool abs_phys;

	vector<string> phys_value_output_file(P_NUM_BUILT_IN_PHYS_VALUE + 1);

#ifdef _mpi_use
	printf("Initializing MPI environment...\n");
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p_num_check);
#else
	my_rank = 0;
	p_num_check = 1;
#endif

	phys_value_output_file[0] = "E";
	phys_value_output_file[1] = "Mx";
	phys_value_output_file[2] = "My";
	phys_value_output_file[3] = "Mz";
	phys_value_output_file[4] = "Q3z2-r2";
	phys_value_output_file[5] = "Qx2-y2";
	phys_value_output_file[6] = "Qxy";
	phys_value_output_file[7] = "Qxz";
	phys_value_output_file[8] = "Qyz";
	phys_value_output_file[9] = "Q2";
	phys_value_output_file[10] = "M2";
	phys_value_output_file[11] = "ss";
	phys_value_output_file[12] = "SG";

	if (my_rank == 0)
		printf("Loading parameters...\n");
	//Loading parameters
	ipp = fopen("param.def", "r");
	fscanf(ipp, "%s %d\n", t_char, &num_pdef);
	fscanf(ipp, "%s %d\n", t_char, &calc_mode);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.load_config);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.dump_config);
	fscanf(ipp, "%s %d\n", t_char, &dump_phys_value);
	fscanf(ipp, "%s %d\n", t_char, &index_start_run);
	fscanf(ipp, "%s %d\n", t_char, &num_run);
	fscanf(ipp, "%s %d\n", t_char, &Ns);
	fscanf(ipp, "%s %d\n", t_char, &num_MC_step);
	fscanf(ipp, "%s %d\n", t_char, &num_MC_step_thermalization);
	fscanf(ipp, "%s %d\n", t_char, &MC_step_per_ave);
	fscanf(ipp, "%s %d\n", t_char, &MC_step_per_ave_ss);
	fscanf(ipp, "%s %lf\n", t_char, &mc_update_param.max_beta_single_spin_flip);
	fscanf(ipp, "%s %lf\n", t_char, &mc_update_param.fluctuation_random_walk);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.num_overrelaxation);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.ex_interval);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.num_inner_process_ex);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.num_ex_opt);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.ex_opt_interval);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.ex_opt_min_beta_index);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.ex_opt_max_beta_index);
	fscanf(ipp, "%s %lf\n", t_char, &mc_update_param.min_beta_loop_move);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.max_num_site_loop);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.num_site_ising_axis);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.loop_flip_method);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.num_cluster_update);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.cluster_update_algorithm);
	fscanf(ipp, "%s %d\n", t_char, &num_replica);
	fscanf(ipp, "%s %d\n", t_char, &num_seperated_replica_set);
	fscanf(ipp, "%s %d\n", t_char, &mc_update_param.auto_correlation);
	fscanf(ipp, "%s %d\n", t_char, &output_phys_value_dynamics);
	fscanf(ipp, "%s %lu\n", t_char, &i_seed);
	fscanf(ipp, "%s %d\n", t_char, &nine_comp_ex_mode);
	fscanf(ipp, "%s %d\n", t_char, &i_abs_phys);
	if(i_abs_phys != 0) {
		abs_phys = true;
	} else {
		abs_phys = false;
	}
	i_clock = time(NULL);
	i_seed = i_clock + my_rank * i_seed;
	dsfmt_init_gen_rand(&dsfmt, (unsigned long) i_seed);
	rtmp = 0.0;
	for (int i = 0; i < 1000000; i++)
	{
		rtmp += dsfmt_genrand_close_open(&dsfmt);
	}
	i_seed += (int) rtmp / 1000000;
	dsfmt_init_gen_rand(&dsfmt, (unsigned long) i_seed);
	printf("my_rank=%d i_clock=%lu\n", my_rank, i_clock);
	printf("my_rank=%d i_seed=%lud\n", my_rank, i_seed);

#ifdef _mpi_use
	if(p_num_check%num_pdef != 0){
		cout << "p_num_check%num_pdef != 0." << endl;
		MPI_Finalize();
		exit(1);
	}
	my_world_size = p_num_check/num_pdef;
	my_world = my_rank/my_world_size;
	if(num_pdef > 1){
		if(my_rank==0)
			cout << "Splitting MPI_COMM_WORLD into " << num_pdef << " parts." << endl;
		MPI_Comm_split(MPI_COMM_WORLD, my_world, my_rank, &MY_MPI_COMM_WORLD);
		MPI_Comm_rank(MY_MPI_COMM_WORLD, &my_splitted_rank);
	}else{
		MY_MPI_COMM_WORLD = MPI_COMM_WORLD;
		my_splitted_rank = my_rank;
		my_world = 0;
		my_world_size = p_num_check;
	}
#else
	my_splitted_rank = my_rank;
	my_world = 0;
	my_world_size = p_num_check;
#endif

	//Loading Hamiltonian
	{
		std::ostringstream oss;
		oss << "two-body-int" << my_world << ".def";
		FILE *Jfp, *Ifp, *Sfp, *Mfp;

		if(my_splitted_rank == 0)
			cout << "Opening... " << oss.str() << endl;
		Jfp = fopen(oss.str().c_str(), "r");
                if (Jfp==NULL) {std::cerr<<oss.str()<<".def cannot be opened!"<<std::endl;exit(1);}

		if(my_splitted_rank == 0)
			cout << "Opening... ia.def..." << endl;
		Ifp = fopen("ia.def", "r");
                if (Ifp==NULL) {std::cerr<<"ia.def cannot be opened!"<<std::endl;exit(1);}

		if(my_splitted_rank == 0)
			cout << "Opening... spin.def..." << endl;
		Sfp = fopen("spin.def", "r");
                if (Sfp==NULL) {std::cerr<<"spin.def cannot be opened!"<<std::endl;exit(1);}

		if(my_splitted_rank == 0)
			cout << "Opening... mag.def..." << endl;
		Mfp = fopen("mag.def", "r");
                if (Mfp==NULL) {std::cerr<<"mag.def cannot be opened!"<<std::endl;exit(1);}


		if(my_splitted_rank == 0)
			cout << "Loading..." << endl;
                if(nine_comp_ex_mode!=0){
			pham = new Ham(Ns, Jfp, Ifp, Sfp, Mfp, true);
		}else{
			pham = new Ham(Ns, Jfp, Ifp, Sfp,  Mfp,false);
		}

		fclose(Jfp);
		fclose(Ifp);
		fclose(Sfp);
	}

	//Loading temperatures
	if (my_rank == 0)
		printf("Loading beta.def...\n");
	iJp = fopen("beta.def", "r");
	fscanf(iJp, "%d", &num_beta);
	num_beta_per_process = num_beta / my_world_size;
	if( !(num_beta_per_process >=1 && num_beta%my_world_size ==0) ){
#ifdef _mpi_use
		MPI_Finalize();
#endif
		exit(1);
	}
	d_malloc1(beta, num_beta);
	for (int i = 0; i <= num_beta - 1; i++)
	{
		fscanf(iJp, "%lf", &beta[i]);
	}
	fclose(iJp);

	//Loading temperatures for reweighting
	if (my_rank == 0)
		printf("Loading beta-rw.def...\n");
	iJp = fopen("beta-rw.def", "r");
	fscanf(iJp, "%d", &num_beta_rw);
	d_malloc1(beta_rw, num_beta_rw);
	for (int i = 0; i <= num_beta_rw - 1; i++)
	{
		fscanf(iJp, "%lf", &beta_rw[i]);
	}
	fclose(iJp);

	//Loading list of user-defined sublattice magnetization
	if (my_rank == 0)
		cout << "Loading sbl-mag.def..." << endl;
	iJp = fopen("sbl-mag.def", "r");
	fscanf(iJp, "%d", &num_sbl_mag);
	sbl_mag_coeff.resize(Ns*num_sbl_mag);
	for (int i=0; i<=Ns*num_sbl_mag-1; ++i){
		double re, im;
		fscanf(iJp, "%lf %lf", &re, &im);
		sbl_mag_coeff[i] = complex<double>(re, im);
	}
	fclose(iJp);

        //Loading cross-product order parameter
        {
	        if (my_rank==0) cout << "Loading cross-product.def..." << endl;
                ifstream input_s("cross-product.def");
                if(!input_s.is_open()) {
                    std::cout<<"cross-product.def is not open or not found! "<<std::endl;
                } else {
                    int num_cross_product;
                    input_s >> num_cross_product;
                    cross_product_coeff.resize(num_cross_product);
                    for (int i=0; i<num_cross_product; ++i) {
                        int num_link;
                        input_s >> num_link;
                        cross_product_coeff[i].resize(num_link);
                        for (int il=0; il<num_link; ++il) {
                            int site1, site2;
                            double re;
                            input_s >> site1;
                            input_s >> site2;
                            input_s >> re;
                            cross_product_coeff[i][il] = boost::make_tuple(site1,site2,re);
                        }
                    }
                }
        }

	//Loading list of spin-spin correlations
	if (my_rank == 0)
		printf("Loading spin-spin-correlation.def...\n");
	iJp = fopen("spin-spin-correlation.def", "r");
	fscanf(iJp, "%d", &num_ss);
	i_malloc2(ss_pair, num_ss, 2)
	for (int i = 0; i <= num_ss - 1; i++)
	{
		fscanf(iJp, "%d %d", &(ss_pair[i][0]), &(ss_pair[i][1]));
	}
	fclose(iJp);

	if (my_rank == 0)
	{
		printf("num_run=%i\n", num_run);
		printf("Ns=%i\n", pham->Ns);
		printf("num_two_body_int=%i\n", pham->num_two_body_int);
		printf("num_beta=%i\n", num_beta);
		printf("num_beta_rw=%i\n", num_beta_rw);
		printf("num_ia=%i\n", pham->num_ia);
		printf("num_process=%i\n", my_world_size);
		printf("fluctuation_random_walk=%lf\n", mc_update_param.fluctuation_random_walk);
		printf("num_beta_per_process=%i\n", num_beta_per_process);
		printf("num_replica=%i\n", num_replica);
		printf("num_seperated_replica_set=%i\n", num_seperated_replica_set);
		printf("num_cluster_update=%d\n", mc_update_param.num_cluster_update);
		printf(
				"num_MC_Step=%d, num_MC_step_thermalization=%d, MC_step_per_ave=%d\n",
				num_MC_step, num_MC_step_thermalization, MC_step_per_ave);
		printf("num_ss=%d\n", num_ss);
		printf("num_sbl_mag=%d\n", num_sbl_mag);
	}

	//Initial MC 
	if(my_splitted_rank == 0)
	{
		pvr.resize(num_beta);
		sr.resize(num_beta);
		ssr.resize(num_beta);
	}
	t_pvr.resize(num_beta_per_process);
	t_sr.resize(num_beta_per_process);
	t_ssr.resize(num_beta_per_process);

	d_malloc2(ex_rate, num_beta_per_process, 10)
	d_malloc1(t_q2_dynamics, num_MC_step);

	utime.resize(5);

	for(int i_run = index_start_run; i_run < index_start_run + num_run; ++ i_run){
		cout << "Staring i_run=" << i_run << endl;
		if(calc_mode == 0){
			printf("my_rank=%d Calling mc...\n", my_rank);
			try{
				int status = mc(my_world, i_run, &dsfmt, pham, num_beta, beta,
                                        num_sbl_mag, sbl_mag_coeff, 
                                        cross_product_coeff, 
					num_replica, num_seperated_replica_set,
					num_MC_step, num_MC_step_thermalization, MC_step_per_ave, MC_step_per_ave_ss,
					&mc_update_param, ex_rate, t_q2_dynamics, t_pvr, t_sr, num_ss, ss_pair, t_ssr, auto_correlation, utime,abs_phys);
				if(status == -1){
					cout << "Error occurs in mc" << endl;
#ifdef _mpi_use
					MPI_Finalize();
#endif
					exit(1);
				}
			}
			catch(const char* str){
				cout << str;
				exit(1);
			}
			printf("my_rank=%d Exiting mc...\n", my_rank);

			/**** gather physical quantities from all processes ****/
			printf("my_rank=%d Gathering phys...\n", my_rank);
			#ifdef _mpi_use
			printf("  my_rank=%d Gathering default phys...\n", my_rank);
			vector<PhysValueRecorder>::iterator input_it = t_pvr.begin();
			vector<PhysValueRecorder>::iterator output_it;
			if(my_splitted_rank == 0)
				output_it = pvr.begin();
			tmpl_mpi_gather_phys_value<vector<PhysValueRecorder>::iterator, PhysValueRecorder>
				(MY_MPI_COMM_WORLD, input_it, output_it, num_beta, my_world_size, my_splitted_rank);

			printf("  my_rank=%d Gathering SG phys...\n", my_rank);
			vector<PhysValueRecorder>::iterator input_it2 = t_sr.begin();
			vector<PhysValueRecorder>::iterator output_it2;
			if(my_splitted_rank == 0)
				output_it2 = sr.begin();
			tmpl_mpi_gather_phys_value<vector<PhysValueRecorder>::iterator, PhysValueRecorder>
				(MY_MPI_COMM_WORLD, input_it2, output_it2, num_beta, my_world_size, my_splitted_rank);

			printf("  my_rank=%d Gathering spin-spin correlations ...\n", my_rank);
			vector<PhysValueRecorder>::iterator input_it3 = t_ssr.begin();
			vector<PhysValueRecorder>::iterator output_it3;
			if(my_splitted_rank == 0)
				output_it3 = ssr.begin();
			tmpl_mpi_gather_phys_value<vector<PhysValueRecorder>::iterator, PhysValueRecorder>
				(MY_MPI_COMM_WORLD, input_it3, output_it3, num_beta, my_world_size, my_splitted_rank);

			#else
			for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
			{
				pvr[i_beta] = t_pvr[i_beta];
				sr[i_beta] = t_sr[i_beta];
				ssr[i_beta] = t_ssr[i_beta];
			}
			#endif

			printf("my_rank=%d Writing log...\n", my_rank);
			{
				std::ostringstream oss;
				oss << "log-run" << i_run << "-" << my_rank;
				iJp = fopen(oss.str().c_str(), "w");
				for (i_beta = 0; i_beta <= num_beta_per_process - 1; i_beta++)
				{
					fprintf(iJp,
							"i_beta= %d EXCHANGE= %25.16lf SINGLE_SPIN= %f OVER_RELAX= %f NUM_SITE_VISITED= %f NUM_TRY= %f NUM_LOOP_FOUND= %f NUM_LOOP_FLIPPED= %f NUM_SITE_FLIPPED= %f OVERLAP_CLUSTER= %f\n",
							i_beta + num_beta_per_process * my_splitted_rank,
								ex_rate[i_beta][0], ex_rate[i_beta][1], ex_rate[i_beta][2],
							ex_rate[i_beta][4], ex_rate[i_beta][5], ex_rate[i_beta][6], ex_rate[i_beta][7], ex_rate[i_beta][8],
							ex_rate[i_beta][9]
						);
					#ifdef _utime
					fprintf(iJp,
						"utime: i_beta= %d ORC= %f SINGLE_SPIN= %f LOOP_MOVE= %f PHYS= %f EX= %f\n", i_beta + num_beta_per_process * my_splitted_rank,
						utime[0], utime[1], utime[2], utime[3], utime[4]);
					#endif
				}
				fclose(iJp);
			}
	
			if(dump_phys_value && my_splitted_rank == 0)
			{
				std::ostringstream oss[2];
				ofstream ofs[2];

				cout << "saving phys values..." << endl;
				oss[0] << "pvr-def" << my_world << "-run" << i_run << ".dat";
				oss[1] << "sr-def" << my_world << "-run" << i_run << ".dat";
				ofs[0].open(oss[0].str().c_str());
				ofs[1].open(oss[1].str().c_str());
				for(i_beta = 0; i_beta < num_beta; ++ i_beta)
				{
					pvr[i_beta].dump(ofs[0]);
					sr[i_beta].dump(ofs[1]);
				}
			}
		}
		else if(calc_mode == 1 && my_splitted_rank == 0)
		{
			try
			{
				std::ostringstream oss[2];
				ifstream ifs[2];

				cout << "loading phys values from data files ..." << endl;
				oss[0] << "pvr-def" << my_world << "-run" << i_run << ".dat";
				oss[1] << "sr-def" << my_world << "-run" << i_run << ".dat";
				ifs[0].open(oss[0].str().c_str(), ios::in);
				ifs[1].open(oss[1].str().c_str(), ios::in);
				for(i_beta = 0; i_beta < num_beta; ++ i_beta)
				{
					pvr[i_beta] = PhysValueRecorder::load(ifs[0]);
					sr[i_beta] = PhysValueRecorder::load(ifs[1]);
				}
			}
			catch(const char* str)
			{
				cout << str;
				exit(1);
			}
		}

		#ifdef _mpi_use
		MPI_Barrier(MY_MPI_COMM_WORLD);
		#endif
		/**** output physical quantities ****/
		if (my_splitted_rank == 0)
		{
			std::ostringstream oss[2];
			oss[0] << "-def" << my_world << "-run" << i_run << ".dat";
			oss[1] << "-rw-def" << my_world << "-run" << i_run << ".dat";
			try
			{
				vector<vector<vector<double> > > dvec3, ave_SG, ave_ssc;
				dvec3.resize(num_beta);
				ave_SG.resize(num_beta);
				ave_ssc.resize(num_beta);
	
				for(i_beta = 0; i_beta < num_beta; ++ i_beta)
				{
					pvr[i_beta].average(dvec3[i_beta]);
					sr[i_beta].average(ave_SG[i_beta]);
					ssr[i_beta].average(ave_ssc[i_beta]);
				}
	
				for(int i_phys_value = 0; i_phys_value < P_NUM_BUILT_IN_PHYS_VALUE; ++ i_phys_value)
				{
					iJp = fopen((phys_value_output_file[i_phys_value] + oss[0].str()).c_str(), "w");
					for(i_beta = 0; i_beta < num_beta; ++ i_beta){
						fprintf(iJp, "%25.16le", beta[i_beta]);
						int num_power = pvr[i_beta].get_num_power();
						for(int i_power=0; i_power<num_power; ++i_power)
							fprintf(iJp, "%25.16le", dvec3[i_beta][i_power][i_phys_value]);
						fprintf(iJp, "\n");
					}
					fclose(iJp);
				}
	
				for(int i_sbl_mag = 0; i_sbl_mag < num_sbl_mag; ++i_sbl_mag){
					std::ostringstream t_oss;
					t_oss << i_sbl_mag;
					iJp = fopen((string("sbl-mag")+t_oss.str()+oss[0].str()).c_str(), "w");
					for(i_beta = 0; i_beta < num_beta; ++ i_beta){
						fprintf(iJp, "%25.16le", beta[i_beta]);
						int num_power = pvr[i_beta].get_num_power();
						for(int i_power=0; i_power<num_power; ++i_power){
							for(int ix=0; ix<4; ++ix)
								fprintf(iJp, "%25.16le", dvec3[i_beta][i_power][P_NUM_BUILT_IN_PHYS_VALUE+i_sbl_mag*4+ix]);
						}
						fprintf(iJp, "\n");
					}
					fclose(iJp);
				}

				for(int i_cr = 0; i_cr < cross_product_coeff.size(); ++i_cr){
					std::ostringstream t_oss;
                                        t_oss << i_cr;
                                        ofstream output_s(("cross-product-"+t_oss.str()+oss[0].str()+".dat").c_str());
					for(int i_beta= 0; i_beta< num_beta; ++ i_beta){
						output_s << beta[i_beta] << " " ;
						int num_power = pvr[0].get_num_power();
						for(int i_power=0; i_power<num_power; ++i_power){
							for(int ix=0; ix<4; ++ix)
								output_s << dvec3[i_beta][i_power][P_NUM_BUILT_IN_PHYS_VALUE+num_sbl_mag*4+i_cr*4+ix] << " ";
						}
                                                output_s << std::endl;
					}
				}
	
				iJp = fopen((phys_value_output_file[P_NUM_BUILT_IN_PHYS_VALUE] + oss[0].str()).c_str(), "w");
				for(i_beta = 0; i_beta < num_beta; ++ i_beta){
					fprintf(iJp, "%25.16le %25.16le %25.16le\n", beta[i_beta], ave_SG[i_beta][0][0], ave_SG[i_beta][1][0]);
				}
				fclose(iJp);

				//output spin-spin correlation
				for(int i_power=0; i_power < 3; ++i_power){
					std::ostringstream t_oss;
					t_oss << "spin-spin-correlation-pow" << i_power;
					iJp = fopen((t_oss.str() + oss[0].str()).c_str(), "w");
					for(int i_beta = 0; i_beta < num_beta; ++ i_beta){
						fprintf(iJp, "%25.16le", beta[i_beta]);
						for(int i_ss_pair = 0; i_ss_pair < num_ss; ++ i_ss_pair)
							fprintf(iJp, " %25.16le", ave_ssc[i_beta][i_power][i_ss_pair]);
						fprintf(iJp, "\n");
					}
					fclose(iJp);
				}
	
				//reweighted version
				int num_sample_disregarded = num_MC_step_thermalization/MC_step_per_ave;
				dvec3.resize(num_beta_rw);
				ave_SG.resize(num_beta_rw);
				for(int i_beta_rw = 0; i_beta_rw < num_beta_rw; ++ i_beta_rw)
				{
					double t_beta_rw = beta_rw[i_beta_rw];
					double rtmp = fabs(beta[0]-t_beta_rw);
					int i_beta_org = 0;
					for(int i_beta = 0; i_beta < num_beta; ++ i_beta)
					{
						if(fabs(beta[i_beta] - t_beta_rw) < rtmp)
						{
							i_beta_org = i_beta;
							rtmp = fabs(beta[i_beta] - t_beta_rw);
						}
					}
					cout << "Reweighting i_beta_rw= " << i_beta_rw << " using i_beta= " << i_beta_org << endl;
					double d_beta = t_beta_rw - beta[i_beta_org];
					pvr[i_beta_org].reweight(num_sample_disregarded, dvec3[i_beta_rw], d_beta);
					sr[i_beta_org].reweight(num_sample_disregarded, ave_SG[i_beta_rw], d_beta);
				}
	
				for(int i_phys_value = 0; i_phys_value < P_NUM_BUILT_IN_PHYS_VALUE; ++ i_phys_value)
				{
					iJp = fopen( (phys_value_output_file[i_phys_value] + oss[1].str()).c_str(), "w");
					for(int i_beta_rw = 0; i_beta_rw < num_beta_rw; ++ i_beta_rw){
						//cout << i_phys_value << " " << i_beta_rw << endl;
						fprintf(iJp, "%25.16le", beta_rw[i_beta_rw]);
						int num_power = pvr[0].get_num_power();
						for(int i_power=0; i_power<num_power; ++i_power){
							fprintf(iJp, "%25.16le", dvec3[i_beta_rw][i_power][i_phys_value]);
						}
						fprintf(iJp, "\n");
						//fprintf(iJp, "%25.16le %25.16le %25.16le %25.16le\n",
							//beta_rw[i_beta_rw], dvec3[i_beta_rw][0][i_phys_value], dvec3[i_beta_rw][1][i_phys_value], dvec3[i_beta_rw][2][i_phys_value]);
					}
					fclose(iJp);
				}

				for(int i_sbl_mag = 0; i_sbl_mag < num_sbl_mag; ++i_sbl_mag){
					std::ostringstream t_oss;
					t_oss << i_sbl_mag;
					iJp = fopen((string("sbl-mag")+t_oss.str()+oss[1].str()).c_str(), "w");
					for(int i_beta_rw = 0; i_beta_rw < num_beta_rw; ++ i_beta_rw){
						fprintf(iJp, "%25.16le", beta_rw[i_beta_rw]);
						int num_power = pvr[0].get_num_power();
						for(int i_power=0; i_power<num_power; ++i_power){
							for(int ix=0; ix<4; ++ix)
								fprintf(iJp, "%25.16le", dvec3[i_beta_rw][i_power][P_NUM_BUILT_IN_PHYS_VALUE+i_sbl_mag*4+ix]);
						}
						fprintf(iJp, "\n");
					}
					fclose(iJp);
				}

				for(int i_cr = 0; i_cr < cross_product_coeff.size(); ++i_cr){
					std::ostringstream t_oss;
                                        t_oss << i_cr;
                                        ofstream output_s(("cross-product-rw-"+t_oss.str()+oss[1].str()+".dat").c_str());
					for(int i_beta_rw = 0; i_beta_rw < num_beta_rw; ++ i_beta_rw){
						output_s << beta_rw[i_beta_rw] << " " ;
						int num_power = pvr[0].get_num_power();
						for(int i_power=0; i_power<num_power; ++i_power){
							for(int ix=0; ix<4; ++ix)
								output_s << dvec3[i_beta_rw][i_power][P_NUM_BUILT_IN_PHYS_VALUE+num_sbl_mag*4+i_cr*4+ix] << " " ;
						}
                                                output_s << std::endl;
					}
				}
	
				iJp = fopen((phys_value_output_file[P_NUM_BUILT_IN_PHYS_VALUE] + oss[1].str()).c_str(), "w");
				for(int i_beta_rw = 0; i_beta_rw < num_beta_rw; ++ i_beta_rw)
				{
					fprintf(iJp, "%25.16le %25.16le %25.16le\n", beta_rw[i_beta_rw], ave_SG[i_beta_rw][0][0], ave_SG[i_beta_rw][1][0]);
				}
				fclose(iJp);
			}
			catch(const char* str)
			{
				cout << str;
				exit(1);
			}
		}

		/**** output dynamics of M^2, Q^2 and q^2, auto correlations ****/
		if(calc_mode == 0)
		{
			try
			{
				std::ostringstream oss;

				if(output_phys_value_dynamics == 1)
				{
					oss << "phys_value_dynamics-def" << my_world << "-" << my_splitted_rank << "-run" << i_run << ".dat";
					//iJp = fopen(oss.str().c_str(), "w");
                                        ofstream output_s(oss.str().c_str());
		
					PhysValueRecorder &ptr_pvr = t_pvr[num_beta_per_process-1];
					PhysValueRecorder &ptr_sr = t_sr[num_beta_per_process-1];
					int num_sample = ptr_pvr.get_num_sample();
					int num_replica = ptr_pvr.get_num_data_set();
					int num_replica_set = ptr_sr.get_num_data_set();
					
					for (int i_sample = 0; i_sample < num_sample; ++ i_sample)
					{
						double M2 = 0.0;
						for(int i_data_set = 0; i_data_set < num_replica; ++ i_data_set)
							M2 += ptr_pvr(i_data_set, 10, i_sample);
		
						double Q2 = 0.0;
						for(int i_data_set = 0; i_data_set < num_replica; ++ i_data_set)
							Q2 += ptr_pvr(i_data_set, 9, i_sample);
		
						double q2 = 0.0;
						for(int i_data_set = 0; i_data_set < num_replica_set; ++ i_data_set)
							q2 += ptr_sr(i_data_set, 0, i_sample);
						
						//fprintf(iJp, "%d %25.16le %25.16le %25.16le\n",
							//i_sample, M2/num_replica, Q2/num_replica, q2/num_replica_set);
						output_s << i_sample  << "     "<< M2/num_replica  << "     "<< Q2/num_replica  << "     "<< q2/num_replica_set << "     ";
                                                {
                                                    for (int i_phys=0; i_phys<ptr_pvr.get_num_phys_value(); ++i_phys) {
                                                        double tmp = 0.0;
						        for(int i_data_set = 0; i_data_set < num_replica; ++ i_data_set) {
							    tmp += ptr_pvr(i_data_set, i_phys, i_sample);
                                                        }
                                                        tmp /= num_replica;
                                                        output_s << " " << tmp << " ";
                                                    }
                                                }
                                                output_s << std::endl;
					}
					//fclose(iJp);
				}

				if(mc_update_param.auto_correlation != 0)
				{
					oss.str("");
					oss << "auto-correlation-def" << my_world << "-" << my_splitted_rank << "-run" << i_run << ".dat";
					iJp = fopen(oss.str().c_str(), "w");
	
					int num_sample = num_MC_step - num_MC_step_thermalization;
					for(int i_sample = 0; i_sample < num_sample; ++ i_sample)
					{
						for(int i_beta = 0; i_beta < num_beta_per_process; ++ i_beta)
							fprintf(iJp, "%25.16le", auto_correlation[i_beta][i_sample]);
						fprintf(iJp, "\n");
					}
					fclose(iJp);
				}
			}
			catch(const char* str)
			{
				cout << str;
				exit(1);
			}
		}
	}

	if(mc_update_param.num_ex_opt > 0 && my_rank == 0)
	{
		iJp = fopen("beta_opt.def", "w");
		fprintf(iJp, "%d\n", num_beta);
		for(int i_beta = 0; i_beta < num_beta; ++ i_beta)
			fprintf(iJp, "%25.16le\n", beta[i_beta]);
		fclose(iJp);
	}

	delete pham;
	d_free1(beta, num_beta)
	d_free2(ex_rate, num_beta_per_process, 2)
	d_free1(t_q2_dynamics, num_MC_step)
	i_free2(ss_pair, num_ss, 2)

#ifdef _mpi_use
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
}
