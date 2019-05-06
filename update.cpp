#include "update.h"

double calc_BE(double S0, double S1, double S2, double B00, double B01,
		double B02, double B11, double B12, double B22)
{
	return -S0 * (B00 * S0 + 2.0 * B01 * S1 + 2.0 * B02 * S2) - S1 * (B11 * S1
			+ 2.0 * B12 * S2) - S2 * B22 * S2;
}

double calc_dE(double S0, double S1, double S2, double new_S0,
		double new_S1, double new_S2, double *local_h, double B00, double B01,
		double B02, double B11, double B12, double B22)
{
	return -local_h[0] * (new_S0 - S0) - local_h[1] * (new_S1 - S1)
			- local_h[2] * (new_S2 - S2) - calc_BE(S0, S1, S2, B00, B01, B02,
			B11, B12, B22) + calc_BE(new_S0, new_S1, new_S2, B00, B01, B02,
			B11, B12, B22);
}

void calc_local_field(int i_site, Ham *ham, double **config,
		int *num_nn_two_body_int, int **idx_cnct_bond,
		int **idx_cnctd_site, double **ia, double *local_h, double *B00,
		double *B01, double *B02, double *B11, double *B12, double *B22)
{
	int i_nn;
	int i_site2, i_coordinate,i_two_body_int;
	double rtmp1;
	bool nine_comp_ex_mode = ham->nine_comp_ex_mode;

	//bi-linear
	local_h[0] = 0.0;
	local_h[1] = 0.0;
	local_h[2] = 0.0;

	//bi-quadratic
	*B00 = 0.0;
	*B01 = 0.0;
	*B02 = 0.0;
	*B11 = 0.0;
	*B12 = 0.0;
	*B22 = 0.0;
	if(nine_comp_ex_mode){
		for (i_nn = 0; i_nn <= num_nn_two_body_int[i_site] - 1; i_nn++)
		{
			i_two_body_int = idx_cnct_bond[i_site][i_nn]%ham->num_two_body_int;
			i_site2 = idx_cnctd_site[i_site][i_nn];

			double J00,J01,J02,J10,J11,J12,J20,J21,J22;

			if(idx_cnct_bond[i_site][i_nn] < ham->num_two_body_int) 
			{
				J00 = ham->two_body_int[i_two_body_int][0];
				J01 = ham->two_body_int[i_two_body_int][1];
				J02 = ham->two_body_int[i_two_body_int][2];
				J10 = ham->two_body_int[i_two_body_int][3];
				J11 = ham->two_body_int[i_two_body_int][4];
				J12 = ham->two_body_int[i_two_body_int][5];
				J20 = ham->two_body_int[i_two_body_int][6];
				J21 = ham->two_body_int[i_two_body_int][7];
				J22 = ham->two_body_int[i_two_body_int][8];
			}else{
				J00 = ham->two_body_int[i_two_body_int][0];
				J10 = ham->two_body_int[i_two_body_int][1];
				J20 = ham->two_body_int[i_two_body_int][2];
				J01 = ham->two_body_int[i_two_body_int][3];
				J11 = ham->two_body_int[i_two_body_int][4];
				J21 = ham->two_body_int[i_two_body_int][5];
				J02 = ham->two_body_int[i_two_body_int][6];
				J12 = ham->two_body_int[i_two_body_int][7];
				J22 = ham->two_body_int[i_two_body_int][8];
			}

			double S0 = config[i_site2][0];
			double S1 = config[i_site2][1];
			double S2 = config[i_site2][2];

			local_h[0] += J00*S0+J01*S1+J02*S2;
			local_h[1] += J10*S0+J11*S1+J12*S2;
			local_h[2] += J20*S0+J21*S1+J22*S2;

			rtmp1 = ham->two_body_int[i_two_body_int][9];
			*B00 += rtmp1 * S0 * S0;
			*B01 += rtmp1 * S0 * S1;
			*B02 += rtmp1 * S0 * S2;
			*B11 += rtmp1 * S1 * S1;
			*B12 += rtmp1 * S1 * S2;
			*B22 += rtmp1 * S2 * S2;
		}
	}else{
		for (i_nn = 0; i_nn <= num_nn_two_body_int[i_site] - 1; i_nn++)
		{
			i_two_body_int = idx_cnct_bond[i_site][i_nn];
			i_site2 = idx_cnctd_site[i_site][i_nn];

			for (i_coordinate = 0; i_coordinate < 3; i_coordinate++)
			{
				local_h[i_coordinate]
						+= ham->two_body_int[i_two_body_int][i_coordinate]
								* config[i_site2][i_coordinate];
			}


			rtmp1 = ham->two_body_int[i_two_body_int][3];
			*B00 += rtmp1 * config[i_site2][0] * config[i_site2][0];
			*B01 += rtmp1 * config[i_site2][0] * config[i_site2][1];
			*B02 += rtmp1 * config[i_site2][0] * config[i_site2][2];
			*B11 += rtmp1 * config[i_site2][1] * config[i_site2][1];
			*B12 += rtmp1 * config[i_site2][1] * config[i_site2][2];
			*B22 += rtmp1 * config[i_site2][2] * config[i_site2][2];
		}
	}

	//Ising anisotropy
	rtmp1 = ia[i_site][0];
	*B00 += rtmp1 * ia[i_site][1] * ia[i_site][1];
	*B01 += rtmp1 * ia[i_site][1] * ia[i_site][2];
	*B02 += rtmp1 * ia[i_site][1] * ia[i_site][3];
	*B11 += rtmp1 * ia[i_site][2] * ia[i_site][2];
	*B12 += rtmp1 * ia[i_site][2] * ia[i_site][3];
	*B22 += rtmp1 * ia[i_site][3] * ia[i_site][3];
}

#ifndef _ising
void update(dsfmt_t *dsfmt, Replica *p_replica, double beta,
		Ham *ham, int **idx_cnct_bond,
		int **pointer_to_counterpart_Jij, int *coord_num,
		double **ia, double **mag, int num_OR, unsigned long *update_counter, double t_fluctuation_random_walk, bool dry_run)
{
	int i_site, i_site2, i, j, i_coordinate, i_counterpart, i_Jij;
	double dE, rtmp1, rtmp2, rtmp3, rtmp4;
	double local_h[3], local_h_normalized[3], B00, B01, B02, B11, B12, B22;
	double new_config[3], inv_norm_local_h;
	double S0, S1, S2;
	const double PI = 3.1415926535;

	for (i_site = 0; i_site <= ham->Ns - 1; i_site++)
	{
		/**** calculate local field ****/
		calc_local_field(i_site, ham, p_replica->config, coord_num,
				idx_cnct_bond, pointer_to_counterpart_Jij, ia, local_h, &B00,
				&B01, &B02, &B11, &B12, &B22);
                for (int i=0; i < 3; ++i)
                        local_h[i] += mag[i_site][i];

		S0 = p_replica->config[i_site][0];
		S1 = p_replica->config[i_site][1];
		S2 = p_replica->config[i_site][2];

		if (num_OR > 0 && ham->spin[i_site] != ISING_SPIN)
		{
			/**** over-relaxation around local field ****/
			inv_norm_local_h = 1.0 / sqrt(norm2(local_h));
			local_h_normalized[0] = local_h[0] * inv_norm_local_h;
			local_h_normalized[1] = local_h[1] * inv_norm_local_h;
			local_h_normalized[2] = local_h[2] * inv_norm_local_h;

			rtmp1 = 2.0 * (local_h_normalized[0] * S0 + local_h_normalized[1]
					* S1 + local_h_normalized[2] * S2);
			new_config[0] = rtmp1 * local_h_normalized[0] - S0;
			new_config[1] = rtmp1 * local_h_normalized[1] - S1;
			new_config[2] = rtmp1 * local_h_normalized[2] - S2;

			dE = calc_dE(S0, S1, S2, new_config[0], new_config[1],
					new_config[2], local_h, B00, B01, B02, B11, B12, B22);
			if (exp(-beta * dE) > dsfmt_genrand_close_open(dsfmt))
			{
				S0 = new_config[0];
				S1 = new_config[1];
				S2 = new_config[2];
				update_counter[1]++;
				p_replica->E += dE;
			}
		}

		/**** propose new configuration (Metropolis) ****/
		if (ham->spin[i_site] == ISING_SPIN) //Ising
		{
			double dsign = 2.0 * (dsfmt_genrand_uint32(dsfmt) % 2) - 1.0;
			new_config[0] = dsign * S0;
			new_config[1] = dsign * S1;
			new_config[2] = dsign * S2;

		}
		else if (ham->spin[i_site] == HEISENBERG_SPIN)
		{ //Heisenberg
                    if (t_fluctuation_random_walk<1.0) {
			rtmp1 = 2.0 * PI * t_fluctuation_random_walk
					* (dsfmt_genrand_close_open(dsfmt) - 0.5); //xy
			rtmp2 = S2 + 2.0 * t_fluctuation_random_walk
					* (dsfmt_genrand_close_open(dsfmt) - 0.5); //z
			rtmp2 = 0.5 * (rtmp2 + 1.0);
			if (rtmp2 >= 0)
			{
				new_config[2] = 2.0 * (rtmp2 - ((int) rtmp2)) - 1.0;
			}
			else
			{
				new_config[2] = 2.0 * (1.0 + rtmp2 - ((int) rtmp2)) - 1.0;
			}
			rtmp2 = sin(rtmp1);
			rtmp3 = cos(rtmp1);
			rtmp4 = sqrt((1 - new_config[2] * new_config[2]) / (S0 * S0 + S1
					* S1));
			new_config[0] = rtmp4 * (rtmp3 * S0 - rtmp2 * S1);
			new_config[1] = rtmp4 * (rtmp2 * S0 + rtmp3 * S1);
                    } else {
			while(1)
			{
				rtmp1 = 1.0-2.0*dsfmt_genrand_open_close(dsfmt);
				rtmp2 = 1.0-2.0*dsfmt_genrand_open_close(dsfmt);
				rtmp3 = rtmp1*rtmp1+rtmp2*rtmp2;
				if(rtmp3<1.0) break;
			}
			rtmp4 = sqrt(1.0-rtmp3);
			new_config[0] = 2.0*rtmp1*rtmp4;
			new_config[1] = 2.0*rtmp2*rtmp4;
			new_config[2] = (1.0-2.0*rtmp3);
			//cout << "i_site= " << " x= " << S0 << " y= " << S1 << " z= " << S2 << endl;
                    }
		}
		else if (ham->spin[i_site] == XY_SPIN)
		{ //XY spins
                        double theta = 2 * M_PI * dsfmt_genrand_open_close(dsfmt);
			new_config[0] = std::cos(theta);
			new_config[1] = std::sin(theta);
			new_config[2] = 0.0;
                }
		else {
			throw std::runtime_error("Unknown spin type");
                }

		/**** calculate energy change ****/
		dE = calc_dE(S0, S1, S2, new_config[0], new_config[1], new_config[2],
				local_h, B00, B01, B02, B11, B12, B22);

		/**** update energy, configurations ****/
		if (exp(-beta * dE) > dsfmt_genrand_close_open(dsfmt))
		{
			update_counter[0]++;
			S0 = new_config[0];
			S1 = new_config[1];
			S2 = new_config[2];
			p_replica->E += dE;
		}

		p_replica->config[i_site][0] = S0;
		p_replica->config[i_site][1] = S1;
		p_replica->config[i_site][2] = S2;
	}
}
#endif

#ifdef _ising
void update(dsfmt_t *dsfmt, Replica *p_replica, double beta,
		Ham *ham, int **idx_cnct_bond,
		int **idx_cnctd_site, int *num_nn_two_body_int,
		double **ia, double **mag, int num_OR, unsigned long *update_counter, bool dry_run)
{
	int i_site, i_site2, i_nn, i_two_body_int;
	double dE;
	double new_config;

	for (i_site = 0; i_site <= ham->Ns - 1; i_site++)
	{
		/**** calculate local field ****/
		new_config = 2.0 * (dsfmt_genrand_uint32(dsfmt) % 2) - 1.0;

		/**** calculate energy change ****/
		dE = 0.0;
		for (i_nn = 0; i_nn <= num_nn_two_body_int[i_site] - 1; i_nn++)
		{
			i_two_body_int = idx_cnct_bond[i_site][i_nn];
			i_site2 = idx_cnctd_site[i_site][i_nn];

			dE -= ham->two_body_int[i_two_body_int][2] * p_replica->config[i_site2][2];
		}
		dE *= new_config - p_replica->config[i_site][2];

		/**** update energy, configurations ****/
		if (exp(-beta * dE) > dsfmt_genrand_close_open(dsfmt))
		{
			update_counter[0]++;
			p_replica->config[i_site][2] = new_config;
			p_replica->E += dE;
		}
	}
}
#endif
