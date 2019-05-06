#include "phys_value.h"

double calc_E(Ham *ham, double **config)
{
	int Ns, num_two_body_int;
	int i, i_coordinate,i_coordinate2;
	double E, **two_body_int;
	int **two_body_int_pair;

	bool nine_comp_ex_mode = ham->nine_comp_ex_mode;

	Ns = ham->Ns;

	E = 0.0;

	num_two_body_int = ham->num_two_body_int;
	two_body_int = ham->two_body_int;
	two_body_int_pair = ham->two_body_int_pair;
	if(nine_comp_ex_mode)
	{
		for (i = 0; i <= num_two_body_int - 1; i++)
		{
			int j=0;
			for (i_coordinate = 0; i_coordinate <= 2; i_coordinate++)
			{
				for (i_coordinate2 = 0; i_coordinate2 <= 2; i_coordinate2++)
				{
					E -= two_body_int[i][j]
						* config[two_body_int_pair[i][0]][i_coordinate]
						* config[two_body_int_pair[i][1]][i_coordinate2];
					j++;
				}
			}

			E -= two_body_int[i][9] * SQUARE(config[two_body_int_pair[i][0]][0]
					* config[two_body_int_pair[i][1]][0]
					+ config[two_body_int_pair[i][0]][1]
							* config[two_body_int_pair[i][1]][1]
					+ config[two_body_int_pair[i][0]][2]
							* config[two_body_int_pair[i][1]][2]);
		}
	}else{
		for (i = 0; i <= num_two_body_int - 1; i++)
		{
			for (i_coordinate = 0; i_coordinate <= 2; i_coordinate++)
			{
				E -= two_body_int[i][i_coordinate]
						* config[two_body_int_pair[i][0]][i_coordinate]
						* config[two_body_int_pair[i][1]][i_coordinate];
			}
			E -= two_body_int[i][3] * SQUARE(config[two_body_int_pair[i][0]][0]
					* config[two_body_int_pair[i][1]][0]
					+ config[two_body_int_pair[i][0]][1]
							* config[two_body_int_pair[i][1]][1]
					+ config[two_body_int_pair[i][0]][2]
							* config[two_body_int_pair[i][1]][2]);
		}
	}

	for (i = 0; i <= ham->num_ia - 1; i++)
	{
		E += ham->ia_value_list[i][0] * (1.0 - SQUARE(
				ham->ia_value_list[i][1] * config[ham->ia_site_list[i]][0]
				+ ham->ia_value_list[i][2] * config[ham->ia_site_list[i]][1]
				+ ham->ia_value_list[i][3] * config[ham->ia_site_list[i]][2]));
	}

	for (int isite = 0; isite < Ns; isite++)
        {
            for (int i = 0; i < 3; ++i)
                E -= ham->mag[isite][i]*config[isite][i];
        }

	return E;
}

/*-------------------------------------------------------------
 * phys_value: [0][0:1]    <|Mx|>, error of <|Mx|>
 *             [1][0:1]    <|My|>, error of <|My|>
 *             [2][0:1]    <|Mz|>, error of <|Mz|>
 *             [3][0:1]    <Q_{3z^2-r-2}>
 *             [4][0:1]    <Q_{x^2-y2}>
 *             [5][0:1]    <Q_{xy}>
 *             [6][0:1]    <Q_{xz}>
 *             [7][0:1]    <Q_{yz}>
 *             [8][0:1]    <Q^2>
 *             [9][0:1]    <M^2>
 *
 * phys_value_suscpt: susceptibility and error of each correspoinding physical quantity
 *-------------------------------------------------------------*/
void calc_phys_value(int Ns, double **config, double *phys_value, int num_ss_pairs, int **ss_pairs,
	int num_sbl_mag, vector<complex<double > > &sbl_mag_coeff, vector<complex<double > > &sbl_mag,
        vector<CROSS_PRODUCT_TYPE>& def_cross_product, vector<complex<double> >& cross_product,  bool abs_phys)
{
	int i, i_site;

	for (i = 0; i <= NUM_PHYS_VALUE - 1; i++)
		phys_value[i] = 0.0;
	for(int i=0; i<3*num_sbl_mag; ++i)
		sbl_mag[i] = 0.0;

	for (i_site = 0; i_site <= Ns - 1; i_site++) {
		for (i = 0; i <= 2; i++)
			phys_value[i] += config[i_site][i];

		phys_value[3] += ((2.0 * SQUARE(config[i_site][2]) - SQUARE( config[i_site][0]) - SQUARE(config[i_site][1])) / sqrt(3.0))/Ns;
		phys_value[4] += (SQUARE(config[i_site][0]) - SQUARE(config[i_site][1]))/Ns;
		phys_value[5] += 2.0 * config[i_site][0] * config[i_site][1]/Ns;
		phys_value[6] += 2.0 * config[i_site][0] * config[i_site][2]/Ns;
		phys_value[7] += 2.0 * config[i_site][1] * config[i_site][2]/Ns;

		for(int i_sbl_mag=0; i_sbl_mag<num_sbl_mag; ++i_sbl_mag){
			for(int ix=0; ix<3; ++ix)
				sbl_mag[ix+3*i_sbl_mag] += sbl_mag_coeff[i_site+i_sbl_mag*Ns]*config[i_site][ix];
		}
	}
	phys_value[8] = SQUARE(phys_value[3]) + SQUARE(phys_value[4]) + SQUARE(
			phys_value[5]) + SQUARE(phys_value[6]) + SQUARE(phys_value[7]);

	if (abs_phys) 
	{
		phys_value[0] = fabs(phys_value[0])/Ns;
		phys_value[1] = fabs(phys_value[1])/Ns;
		phys_value[2] = fabs(phys_value[2])/Ns;
	} else {
		phys_value[0] = (phys_value[0])/Ns;
		phys_value[1] = (phys_value[1])/Ns;
		phys_value[2] = (phys_value[2])/Ns;
	}
	phys_value[9] += SQUARE(phys_value[0]) + SQUARE(phys_value[1]) + SQUARE(phys_value[2]);

	for(int i=0; i<3*num_sbl_mag; ++i) {
		//sbl_mag[i] /= Ns;
		sbl_mag[i] = std::abs(sbl_mag[i])/Ns;
        }

        /*** cross product ***/
        for (int icr=0; icr<def_cross_product.size(); ++icr) {
            std::valarray<double> tmp(0.0, 3);
            std::valarray<double> config1(3), config2(3);
            for (int i_link=0; i_link<def_cross_product[icr].size(); ++i_link) {
                int site1 = (def_cross_product[icr][i_link]).get<0>();
                int site2 = (def_cross_product[icr][i_link]).get<1>();
                double coeff = (def_cross_product[icr][i_link]).get<2>();
                for (int ix=0; ix<3; ++ix) {
                    config1[ix] = config[site1][ix];
                    config2[ix] = config[site2][ix];
                }
                tmp += coeff*cross_prod(config1,config2);
            }
            tmp /= Ns;
            for (int ix=0; ix<3; ++ix) cross_product[3*icr+ix] = std::abs(tmp[ix]);
        }

	{
		double ss = 0.0;
		for(int i_ss_pair=0; i_ss_pair<num_ss_pairs; ++i_ss_pair){
			int i = ss_pairs[i_ss_pair][0];
			int j = ss_pairs[i_ss_pair][1];
			ss += dot_product(config[i], config[j]);
		}
		phys_value[10] = ss;
	}
}

void calc_squared_q(int Ns, Replica *replica1, Replica *replica2,
		double *q2, double *q4)
{
	double rtmp;
	int i_site, i, j;

	double t_q2[9];

	for (i = 0; i < 9; i++)
	{
		t_q2[i] = 0;
	}

	for (i_site = 0; i_site <= Ns - 1; i_site++)
	{
		for (i = 0; i <= 2; i++)
		{
			rtmp = replica1->config[i_site][i];
			for (j = 0; j <= 2; j++)
			{
				t_q2[j + 3 * i] += rtmp * (replica2->config[i_site][j]);
			}
		}
	}

	*q2 = 0.0;
	for (i = 0; i < 9; i++)
	{
		*q2 += t_q2[i] * t_q2[i];
	}
	*q2 /= 1.0 * Ns * Ns;
}

void average(int num_data, double *data, double *ave, double *sigma)
{
	int i_data;
	double rtmp1, rtmp2;

	rtmp1 = 0.0;
	rtmp2 = 0.0;
	for (i_data = 0; i_data <= num_data - 1; i_data++)
	{
		rtmp1 += data[i_data];
		rtmp2 += pow(data[i_data], 2.0);
	};
	*ave = rtmp1 / num_data;
	*sigma = sqrt(rtmp2 / num_data - pow(*ave, 2.0));
}
/****
void calc_spin_correlation(int Ns, double **config, int num_pair, int **pair_list, double *sisj)
{
	for(int i_pair=0; i_pair<num_pair_list; ++i_pair)
	{
		i = pair_list[i_pair][0];
		j = pair_list[i_pair][1];
		sisj[i_pair] = config[i][j];
	}
}
****/
