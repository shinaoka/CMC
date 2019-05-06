/*-------------------------------------------------------------
 * Structure: Hamiltonian
 *-------------------------------------------------------------*/
#ifndef HAM_H
#define HAM_H

#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include "mfmemory.h"
#include "util.h"
#include "spin.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

class Ham {
	public:
		int Ns;

		//Jij
		int num_two_body_int;
		double **two_body_int;
		int **two_body_int_pair;
		vector<bool> flag_nn;
	
		//Ising anisotropy
		int num_ia, *spin;
		double **ia;
		int *ia_site_list;
		double **ia_value_list;

                //magnetic field
		double **mag;

		//for MC update
	 	int *coord_num, **idx_cnct_bond, **idx_cnctd_site;

		bool nine_comp_ex_mode;

		Ham(int Ns, FILE *Jfp, FILE *Ifp, FILE *Sfp, FILE *Mfp, bool nine_comp_ex_mode);
		~Ham();
};
typedef class Ham Ham;

typedef std::vector<boost::tuple<int,int,double> > CROSS_PRODUCT_TYPE;
#endif
