#include "ham.h"

using namespace std;

Ham::Ham(int Ns, FILE *Jfp, FILE *Ifp, FILE *Sfp, FILE *Mfp, bool nine_comp_ex_mode)
{
	int max_len_idx_cnct_bond, i,j;
	//int mfint[7]; //for memory allocator
	std::ostringstream oss;

	this->nine_comp_ex_mode=nine_comp_ex_mode;

	/**** load Hamiltonian ****/
	this->Ns = Ns;
	fscanf(Jfp, "%d", &num_two_body_int);
	fscanf(Ifp, "%d", &num_ia);


	/**** allocate memory  ****/
	//Jij
        if(nine_comp_ex_mode){
		malloc2_seq(two_body_int, num_two_body_int, 10);
	}else{
		malloc2_seq(two_body_int, num_two_body_int, 4);
        }
	malloc2_seq(two_body_int_pair, num_two_body_int, 2);
	flag_nn.resize(num_two_body_int);
	//Ising anistropy
	malloc2_seq(ia, Ns, 4);
	malloc2_seq(ia_value_list, num_ia, 4);
	malloc1_seq(ia_site_list, num_ia);
	malloc1_seq(spin, Ns);
	//magnetic field
	malloc2_seq(mag, Ns, 3);
	//MC updates
	malloc1_seq(coord_num, Ns);

	/**** loading Jij ****/
	for (int i=0; i<=num_two_body_int-1; i++){
		int itmp;
  	if(nine_comp_ex_mode){
			fscanf(Jfp, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &j, 
					&two_body_int_pair[i][0],
					&two_body_int_pair[i][1],
				  &two_body_int[i][0], //Jxx
					&two_body_int[i][1], //Jxy
					&two_body_int[i][2], //Jxz
					&two_body_int[i][3], //Jyx
					&two_body_int[i][4], //Jyy
					&two_body_int[i][5], //Jyz
					&two_body_int[i][6], //Jzx
					&two_body_int[i][7], //Jzy
					&two_body_int[i][8], //Jzz
					&two_body_int[i][9], // biquadratic
					&itmp);
		}else{
			fscanf(Jfp, "%d %d %d %lf %lf %lf %lf %d", &j, &two_body_int_pair[i][0],
					&two_body_int_pair[i][1], &two_body_int[i][0], &two_body_int[i][1],
					&two_body_int[i][2], &two_body_int[i][3], &itmp);
    }
		if (two_body_int_pair[i][0] < 0 || two_body_int_pair[i][0] > Ns - 1)
			cout << "Error in " << oss.str() << endl;
		if (two_body_int_pair[i][1] < 0 || two_body_int_pair[i][1] > Ns - 1)
			cout << "Error in " << oss.str() << endl;
		if(itmp==1){
			flag_nn[i] = true;
		}else if(itmp==0){
			flag_nn[i] = false;
		}else{
			cout << "Invalid flag_nn" << endl;
			exit(1);
		}
	}
	//calculate the coordination number for each site
	for (int i_site = 0; i_site <= Ns-1; i_site++){
		coord_num[i_site] = 0;
	}
	for (int i = 0; i <= num_two_body_int - 1; i++){
		coord_num[two_body_int_pair[i][0]] ++;
		coord_num[two_body_int_pair[i][1]] ++;
	}
	max_len_idx_cnct_bond = 0;
	for (int i_site = 0; i_site <= Ns - 1; i_site++){
		if (max_len_idx_cnct_bond < coord_num[i_site]){
			max_len_idx_cnct_bond = coord_num[i_site];
		}
	}
	malloc2_seq(idx_cnct_bond, Ns, max_len_idx_cnct_bond);
	malloc2_seq(idx_cnctd_site, Ns, max_len_idx_cnct_bond);
	for(int i_site=0; i_site <= Ns-1; i_site++){
		coord_num[i_site] = 0;
		coord_num[i_site] = 0;
	}
	for(int i=0; i<=num_two_body_int-1; i++){
		idx_cnct_bond[two_body_int_pair[i][0]][coord_num[two_body_int_pair[i][0]]] = i;
  	if(nine_comp_ex_mode){
			idx_cnct_bond[two_body_int_pair[i][1]][coord_num[two_body_int_pair[i][1]]] = i+num_two_body_int;
		}else{
			idx_cnct_bond[two_body_int_pair[i][1]][coord_num[two_body_int_pair[i][1]]] = i;
		}
		idx_cnctd_site[two_body_int_pair[i][0]][coord_num[two_body_int_pair[i][0]]]
				= two_body_int_pair[i][1];
		idx_cnctd_site[two_body_int_pair[i][1]][coord_num[two_body_int_pair[i][1]]]
				= two_body_int_pair[i][0];
		coord_num[two_body_int_pair[i][0]] ++;
		coord_num[two_body_int_pair[i][1]] ++;
	}

	/**** load Ising anisotropy ****/
	for(int i = 0; i <= num_ia - 1; i++){
		fscanf(Ifp, "%d %lf %lf %lf %lf", &(ia_site_list[i]),
				&(ia_value_list[i][0]), &(ia_value_list[i][1]), &(ia_value_list[i][2]), &(ia_value_list[i][3]));
		if (ia_site_list[i] < 0 || ia_site_list[i] > Ns - 1)
			printf("Error in ia.def\n");
	}
	/**** magnetic field ****/
	for(int isite = 0; isite <= Ns - 1; isite++){
                int isite2;
		fscanf(Mfp, "%d %lf %lf %lf", &isite2, &(mag[isite][0]), &(mag[isite][1]), &(mag[isite][2]));
                cout << isite << " " << mag[isite][2] << endl;
                if (isite != isite2) {
	        	printf("Error in mag.def\n");
                        exit(1);
                }
        }
	//default easy axis is z axis
	for (int i_site=0; i_site<=Ns-1; i_site++){
		ia[i_site][0] = 0.0;

		ia[i_site][1] = 0.0;
		ia[i_site][2] = 0.0;
		ia[i_site][3] = 1.0;
	}
	//overwrite default easy axis
	for (int i=0; i<=num_ia-1; i++){
		ia[ia_site_list[i]][0] = ia_value_list[i][0]; //alpha
		ia[ia_site_list[i]][1] = ia_value_list[i][1]; //r_x
		ia[ia_site_list[i]][2] = ia_value_list[i][2]; //r_y
		ia[ia_site_list[i]][3] = ia_value_list[i][3]; //r_z
	}
	//normalize
	for (int i_site = 0; i_site <= Ns - 1; i_site++){
		if (fabs(norm2(&ia[i_site][1])) < 1e-8) {
                   throw std::runtime_error("Ising anisotropy axis is invalid!");
                }
		double rtmp = 1.0/sqrt(norm2(&ia[i_site][1]));
		for(int mu = 0; mu < 3; ++ mu)
			ia[i_site][mu+1] *= rtmp;
	}

	for (int i = 0; i <= Ns - 1; i++){
		fscanf(Sfp, "%d", &(spin[i]));
		if (spin[i] != ISING_SPIN && spin[i] != HEISENBERG_SPIN && spin[i] != XY_SPIN)
			printf("Error in spin.def\n");
	}
}

Ham::~Ham(){
	//Jij
	free2_seq(two_body_int);
	free2_seq(two_body_int_pair);
	//Ising anistropy
	free1_seq(spin);
	free2_seq(ia);
	free1_seq(ia_site_list);
	free2_seq(ia_value_list);
	//magnetic field
	free2_seq(mag);
	//MC updates
	free1_seq(coord_num);
	free2_seq(idx_cnct_bond);
	free2_seq(idx_cnctd_site);
}
