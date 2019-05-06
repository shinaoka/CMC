#include "loop_move.h"

int loop_flip(dsfmt_t *dsfmt, int imethod, double **easy_axis_site, int num_update_fix_loop,
	double beta, Replica& replica, int *site_loop, int len_loop,
	Ham *ham, int **idx_cnct_bond, int **idx_cnctd_site, int *coord_num, vector<double> &rwork)
{
	double **config = replica.config;
	//local field
	double local_h[3], B00, B01, B02, B11, B12, B22;
	double dE_tot = 0.0;

	double E[2];

	if(rwork.size() < len_loop*3){
		cout << "Too short rwork!" << endl;
		exit(1);
	}
	//backup the present spin configuration
	for(int idx_loop=0; idx_loop<len_loop; ++idx_loop){
		int i_site = site_loop[idx_loop];
		for(int ix=0; ix<3; ++ix)
			rwork[3*idx_loop+ix] = config[i_site][ix];
	}

	for(int idx_loop=0; idx_loop<len_loop; ++idx_loop){
		int i_site = site_loop[idx_loop];
		double new_config[3];

		//generate new configuration
		if(imethod==0){
			for(int ix=0; ix<3; ++ix)
				new_config[ix] = -config[i_site][ix];
		}else if(imethod==1){
			for(int ix=0; ix<3; ++ix)
				new_config[ix] = rwork[3*mymod(idx_loop+1, len_loop)+ix];
		}else if(imethod==2){
			double *ising_axis = easy_axis_site[i_site];
			double cos = dot_product(ising_axis, config[i_site]);
			for(int ix=0; ix<3; ++ix)
				new_config[ix] = config[i_site][ix] -2.0*cos*ising_axis[ix];
		}

		calc_local_field(i_site, ham, replica.config, coord_num,
			idx_cnct_bond, idx_cnctd_site, ham->ia, local_h, &B00,
			&B01, &B02, &B11, &B12, &B22);
                for (int i=0; i < 3; ++i)
                        local_h[i] += ham->mag[i_site][i];
		dE_tot += calc_dE(
			config[i_site][0], config[i_site][1], config[i_site][2],
			new_config[0], new_config[1], new_config[2],
			local_h, B00, B01, B02, B11, B12, B22);

		//update
		for(int ix=0; ix<3; ++ix)
			config[i_site][ix] = new_config[ix];
	}

	if(dsfmt_genrand_close_open(dsfmt) < exp(-beta*dE_tot)){
		//accept update
		replica.E += dE_tot;
		return len_loop;
	}else{
		//reject update and restore spin configuration
		for(int idx_loop=0; idx_loop<len_loop; ++idx_loop){
			int i_site = site_loop[idx_loop];
			for(int ix=0; ix<3; ++ix)
				config[i_site][ix] = rwork[3*idx_loop+ix];
		}
		return 0;
	}
}

inline void move_to_next_graph(int *i_graph, int *i_nn,
	int ***index_nn_graph, int *num_graph_visited, int *visited_graph, int *graph_chain,
	int *num_site_chain, int *site_chain,
	int *i_site_last_passed)
{
		int i_graph_last_visited;
		int i_graph_next;

		i_graph_next = index_nn_graph[*i_graph][*i_nn][0];

		//move to next graph and register it
		i_graph_last_visited = *i_graph;
		*i_graph = index_nn_graph[i_graph_last_visited][*i_nn][0];
		visited_graph[*i_graph] = *num_graph_visited;
		graph_chain[*num_graph_visited] = *i_graph;
		(*num_graph_visited) ++;

		//adjoint edge site
		*i_site_last_passed = index_nn_graph[i_graph_last_visited][*i_nn][1];
		site_chain[*num_site_chain] = *i_site_last_passed;
		(*num_site_chain) ++;
}

void construct_complete_graph_network(GraphNetwork *graph_n, dsfmt_t *dsfmt, Ham *ham, int **idx_cnct_bond,
			int **idx_cnctd_site, int *coord_num)
{
	int i_bond, i_site, i_site0, i_site1, itmp, **graph_site;
	int i_nn0, i_nn1, i_graph, i, j, *ia_tmp;
	int mfint[7];
	int i_graph0, i_graph1;

	i_malloc2(graph_site, ham->Ns, 3)
	for(i_site = 0; i_site <= ham->Ns - 1; i_site++)
	{
		graph_site[i_site][0] = -1;
		graph_site[i_site][1] = -1;
		graph_site[i_site][2] = 0;
	}
	graph_n->num_graph = 0;
	for(i_bond = 0; i_bond <= ham->num_two_body_int-1; i_bond++)
	{
		if(!ham->flag_nn[i_bond])
			continue;

		i_site0 = ham->two_body_int_pair[i_bond][0];
		i_site1 = ham->two_body_int_pair[i_bond][1];

		itmp = 0;
		for(i = 0; i <= graph_site[i_site0][2] - 1; i++)
		{
			for(j = 0; j <= graph_site[i_site1][2] - 1; j++)
			{
				if(graph_site[i_site0][i] == graph_site[i_site1][j])
				{
					itmp = 1;
					break;
				}
			}
		}
		if(itmp)
			continue;

		graph_site[i_site0][graph_site[i_site0][2]] = graph_n->num_graph;
		graph_site[i_site0][2] ++;
		graph_site[i_site1][graph_site[i_site1][2]] = graph_n->num_graph;
		graph_site[i_site1][2] ++;
		for(i_nn0 = 0; i_nn0 <= coord_num[i_site0]-1; i_nn0++)
		{
			if(!ham->flag_nn[idx_cnct_bond[i_site0][i_nn0]%ham->num_two_body_int])
				continue;

			itmp = idx_cnctd_site[i_site0][i_nn0];
			for(i_nn1 = 0; i_nn1 <= coord_num[i_site1]-1; i_nn1++)
			{
				if(!ham->flag_nn[idx_cnct_bond[i_site1][i_nn1]%ham->num_two_body_int])
					continue;
				if(itmp == idx_cnctd_site[i_site1][i_nn1])
				{
					graph_site[itmp][graph_site[itmp][2]] = graph_n->num_graph;
					graph_site[itmp][2] ++;
				}
			}
		}
		graph_n->num_graph ++;
	}

	/**** count # of sites and nn graph for each graph ****/
	i_calloc1(graph_n->num_site_graph, graph_n->num_graph)
	i_calloc1(graph_n->num_nn_graph, graph_n->num_graph)
	for(i_site = 0; i_site <= ham->Ns - 1; i_site++)
	{
		for(i = 0; i <= graph_site[i_site][2] - 1; i++)
		{
			i_graph = graph_site[i_site][i];
			graph_n->num_site_graph[i_graph] ++;
		}
		if(graph_site[i_site][2] == 2)
		{
			graph_n->num_nn_graph[graph_site[i_site][0]] ++;
			graph_n->num_nn_graph[graph_site[i_site][1]] ++;
		}
	}

	/**** count maximum # of sites and nn graph for each graph ****/
	graph_n->max_num_site_graph = 0;
	graph_n->max_num_nn_graph = 0;
	for(i_graph = 0; i_graph <= graph_n->num_graph - 1; i_graph++)
	{
		if(graph_n->max_num_site_graph < graph_n->num_site_graph[i_graph])
			graph_n->max_num_site_graph = graph_n->num_site_graph[i_graph];
		if(graph_n->max_num_nn_graph < graph_n->num_nn_graph[i_graph])
			graph_n->max_num_nn_graph = graph_n->num_nn_graph[i_graph];
	}

	/**** register sites and nn graph for each graph ****/
	i_malloc2(graph_n->index_site_graph, graph_n->num_graph, graph_n->max_num_site_graph)
	i_malloc3(graph_n->index_nn_graph, graph_n->num_graph, graph_n->max_num_nn_graph, 2)
	for(i_graph = 0; i_graph <= graph_n->num_graph - 1; i_graph++)
	{
		graph_n->num_site_graph[i_graph] = 0;
		graph_n->num_nn_graph[i_graph] = 0;
	}
	for(i_site = 0; i_site <= ham->Ns - 1; i_site++)
	{
		//register site
		for(i = 0; i <= graph_site[i_site][2] - 1; i++)
		{
			i_graph = graph_site[i_site][i];
			graph_n->index_site_graph[i_graph][graph_n->num_site_graph[i_graph]] = i_site;
			graph_n->num_site_graph[i_graph] ++;
		}

		//register nn graph
		if(graph_site[i_site][2] == 2)
		{
			i_graph0 = graph_site[i_site][0];
			i_graph1 = graph_site[i_site][1];

			graph_n->index_nn_graph[i_graph0][graph_n->num_nn_graph[i_graph0]][0] = i_graph1;
			graph_n->index_nn_graph[i_graph1][graph_n->num_nn_graph[i_graph1]][0] = i_graph0;

			graph_n->index_nn_graph[i_graph0][graph_n->num_nn_graph[i_graph0]][1] = i_site;
			graph_n->index_nn_graph[i_graph1][graph_n->num_nn_graph[i_graph1]][1] = i_site;

			graph_n->num_nn_graph[i_graph0] ++;
			graph_n->num_nn_graph[i_graph1] ++;
		}
	}

	/**** count bond in graph ****/
	itmp =  graph_n->max_num_site_graph;
	i_malloc2(graph_n->bond_graph, graph_n->num_graph, itmp*(itmp-1)/2)
	i_malloc1(ia_tmp,graph_n->num_graph);
	for(int ig=0; ig<graph_n->num_graph; ++ig)
		ia_tmp[ig] = 0;

	for(i_bond = 0; i_bond <= ham->num_two_body_int-1; i_bond++)
	{
		itmp = 0;
		i_site0 = ham->two_body_int_pair[i_bond][0];
		i_site1 = ham->two_body_int_pair[i_bond][1];
		for(i = 0; i <= graph_site[i_site0][2] - 1; i ++)
		{
			i_graph0 = graph_site[i_site0][i];
			for(j = 0; j <= graph_site[i_site1][2] - 1; j ++)
			{
				i_graph1 = graph_site[i_site1][j];
				if(i_graph0 == i_graph1)
				{
					graph_n->bond_graph[i_graph0][ia_tmp[i_graph0]] = i_bond;
					ia_tmp[i_graph0] ++;
				}
			}
		}
	}
	i_free1(ia_tmp,graph_n->num_graph);

	i_free2(graph_site, ham->Ns, 3)
}

void free_complete_graph_network(GraphNetwork *graph_n)
{
	int mfint[7], itmp;

	i_free1(graph_n->num_nn_graph, graph_n->num_graph)
	i_free1(graph_n->num_site_graph, graph_n->num_graph)
	i_free2(graph_n->index_site_graph, graph_n->num_graph, graph_n->max_num_site_graph)
	i_free3(graph_n->index_nn_graph, graph_n->num_graph, graph_n->max_num_nn_graph, 2)
	itmp =  graph_n->max_num_site_graph;
	i_free2(graph_n->bond_graph, graph_n->num_graph, itmp*(itmp-1)/2)
}

void allocate_loop_move_work_arrays(LoopMoveWorkArrays *wrk, int Ns, int num_graph, int max_num_nn_graph)
{
	int mfint[7];

	//work arrays (size == ham->Ns)
	i_calloc1(wrk->visited_graph, Ns)
	i_calloc1(wrk->site_chain, Ns)
	i_calloc1(wrk->site_mask, Ns)

	malloc2_seq(wrk->easy_axis_site, Ns, 3);

	//work arrays (size == num_graph)
	i_calloc1(wrk->is_defect_free, num_graph)
	i_calloc1(wrk->index_graph_defect_free, num_graph)
	i_calloc1(wrk->index_defect_free_graph, num_graph)
	i_calloc1(wrk->graph_chain, num_graph)
	//work arrays (size == max_num_nn_graph)
	i_calloc1(wrk->index_nn_graph_visited, max_num_nn_graph)
	//i_calloc1(wrk->index_nn_graph_correctable_defect, max_num_nn_graph)
	i_calloc1(wrk->index_nn_graph_defect_free, max_num_nn_graph)
}

void free_loop_move_work_arrays(LoopMoveWorkArrays *wrk, int Ns, int num_graph, int max_num_nn_graph)
{
	int mfint[7];

	//work arrays (size == ham->Ns)
	i_free1(wrk->visited_graph, Ns)
	i_free1(wrk->site_chain, Ns)
	i_free1(wrk->site_mask, Ns)
	free2_seq(wrk->easy_axis_site);
	//work arrays (size == num_graph)
	i_free1(wrk->is_defect_free, num_graph)
	i_free1(wrk->index_graph_defect_free, num_graph)
	i_free1(wrk->index_defect_free_graph, num_graph)
	i_free1(wrk->graph_chain, num_graph)
	//work arrays (size == max_num_nn_graph)
	i_free1(wrk->index_nn_graph_visited, max_num_nn_graph)
	//i_free1(wrk->index_nn_graph_correctable_defect, max_num_nn_graph)
	i_free1(wrk->index_nn_graph_defect_free, max_num_nn_graph)
}

void loop_move_ISING(dsfmt_t *dsfmt, Replica *p_replica, double beta,
		Ham *ham, int **idx_cnct_bond, int **idx_cnctd_site,
		int *coord_num, double **ia,
		int max_num_site_found, int num_update_fix_loop, int num_site_ising_axis, int loop_flip_method,
		GraphNetwork *graph_n, LoopMoveWorkArrays *wrk,
		int *num_site_visited, int *num_try, int *num_loop_found, int *num_loop_flipped, int *num_site_flipped)
{
	int i_site, i_graph, i_nn, i_flip_start_site_chain, i_flip_end_site_chain, i_affected_graph_start;
	int i, j, k, itmp;
	int i_next_graph, i_next_edge_site, i_loop;
	int i_flip_start_graph_chain, i_flip_end_graph_chain;
	int exit_mode;
	double rtmp, rtmp2;
	const int DEFECT_FOUND=0, VISITED_FOUND=1, FAILED=2;
	double m_site_last_passed, m_next_edge_site;

	int num_nn_graph_visited, num_nn_graph_defect_free, i_site_last_passed;
	int num_graph_defect_free, num_graph_visited, num_site_chain;
	int i_nn_graph_visited_shortest_loop, i_graph_chain_visited_shortest_loop;

	//Ising axis
	double ising_axis[3], ra_tmp[3];
	int i_graph_ising_axis, num_ite_ising_axis = 6;

	vector<double> rwork;

	if(graph_n->num_graph==0)
		return;

	int update_mode;

	int Ns = ham->Ns;

	*num_site_visited = 0;
	*num_loop_found = 0;
	*num_loop_flipped = 0;
	*num_site_flipped = 0;

	rwork.resize(3*Ns);

	double **p_config = p_replica->config;

	/**** estimate easy axis ****/
	set_val_i1(wrk->site_mask, 0, ham->Ns);
	if(num_site_ising_axis == 0)
	{
		for(int i_site = 0; i_site < Ns; ++ i_site)
		{
			wrk->easy_axis_site[i_site][0] = ia[i_site][1];
			wrk->easy_axis_site[i_site][1] = ia[i_site][2];
			wrk->easy_axis_site[i_site][2] = ia[i_site][3];
		}
	}
	else
	{
		for(j = 0; j < num_site_ising_axis; j++)
		{
			wrk->graph_chain[j] = dsfmt_genrand_uint32(dsfmt)%graph_n->num_graph;
		}
		ising_axis[0] = 0.0;
		ising_axis[1] = 0.0;
		ising_axis[2] = 1.0;
		for(k = 0; k < num_ite_ising_axis; k++)
		{
			ra_tmp[0] = 0.0;
			ra_tmp[1] = 0.0;
			ra_tmp[2] = 0.0;
			for(j = 0; j < num_site_ising_axis; j++)
			{
				i_graph = wrk->graph_chain[j];
				for(i = 0; i <= graph_n->num_site_graph[i_graph] - 1; i++)
				{
					i_site = graph_n->index_site_graph[i_graph][i];
					wrk->site_mask[i_site] = 1;
		
					if(dot_product(ising_axis, p_replica->config[i_site]) >= 0.0)
					{
						rtmp = 1.0;
					}
					else
					{
						rtmp = -1.0;
					}
					ra_tmp[0] += rtmp*p_replica->config[i_site][0];
					ra_tmp[1] += rtmp*p_replica->config[i_site][1];
					ra_tmp[2] += rtmp*p_replica->config[i_site][2];
				}
			}
			rtmp = norm2(ra_tmp);
			if(rtmp > 0.0)
			{
				rtmp = 1.0/sqrt(rtmp);
				ising_axis[0] = rtmp*ra_tmp[0];
				ising_axis[1] = rtmp*ra_tmp[1];
				ising_axis[2] = rtmp*ra_tmp[2];
			}
			else
			{
				ising_axis[0] = 0.0;
				ising_axis[1] = 0.0;
				ising_axis[2] = 1.0;
			}
		}
		for(int i_site = 0; i_site < Ns; ++ i_site)
		{
			wrk->easy_axis_site[i_site][0] = ising_axis[0];
			wrk->easy_axis_site[i_site][1] = ising_axis[1];
			wrk->easy_axis_site[i_site][2] = ising_axis[2];
		}
	}

	/**** clean up work arrays and calc magnetization on each graph and site ****/
	num_graph_defect_free = 0;
	for(i_graph = 0; i_graph <= graph_n->num_graph-1; i_graph++)
	{
		wrk->visited_graph[i_graph] = -1;

		int num_site_researved = 0;
		int flux = 0;
		for(i = 0; i < graph_n->num_site_graph[i_graph]; ++ i)
		{
			i_site = graph_n->index_site_graph[i_graph][i];
			num_site_researved += wrk->site_mask[i_site];

			flux += isign(dot_product(p_config[i_site], wrk->easy_axis_site[i_site]));
		}

		if(num_site_researved < graph_n->num_site_graph[i_graph] && abs(flux) == 0)
		{
			wrk->index_graph_defect_free[num_graph_defect_free] = i_graph;
			wrk->index_defect_free_graph[i_graph] = num_graph_defect_free;
			wrk->is_defect_free[i_graph] = 1;
			num_graph_defect_free ++;
		}
		else
		{
			wrk->is_defect_free[i_graph] = 0;
		}
	}
	if(num_graph_defect_free == 0)
		return;

	/**** loop move ****/
	i_loop = 0;
	while(1)
	{
		if(*num_site_visited > max_num_site_found || i_loop > max_num_site_found)
			break;

		num_graph_visited = 0;
		num_site_chain = 0;
		i_graph = wrk->index_graph_defect_free[dsfmt_genrand_uint32(dsfmt)%num_graph_defect_free];
		wrk->visited_graph[i_graph] = num_graph_visited;
		wrk->graph_chain[0] = i_graph;
		i_site_last_passed = -1;
		++ num_graph_visited;

		//search loop
		while(1)
		{
			//search nn site
			num_nn_graph_visited = 0;
			i_nn_graph_visited_shortest_loop = 0;
			i_graph_chain_visited_shortest_loop = -1;
			num_nn_graph_defect_free = 0;
			if(i_site_last_passed != -1)
				m_site_last_passed = dot_product(p_config[i_site_last_passed], wrk->easy_axis_site[i_site_last_passed]);

			for(i_nn = 0; i_nn <= graph_n->num_nn_graph[i_graph]-1; i_nn++)
			{
				i_next_graph = graph_n->index_nn_graph[i_graph][i_nn][0];
				i_next_edge_site = graph_n->index_nn_graph[i_graph][i_nn][1];
				m_next_edge_site = dot_product(p_config[i_next_edge_site], wrk->easy_axis_site[i_next_edge_site]);

				if(i_site_last_passed != -1)
				{
					if(i_next_edge_site == i_site_last_passed)
						continue;
					if(m_site_last_passed*m_next_edge_site > 0.0)
						continue;
				}
	
				if(wrk->site_mask[i_next_edge_site] == 1)
					continue;

				if(wrk->visited_graph[i_next_graph] >= 0 &&
					wrk->visited_graph[i_graph]%2 == 1 - wrk->visited_graph[i_next_graph]%2 )
				{
					wrk->index_nn_graph_visited[num_nn_graph_visited] = i_nn;

					i = wrk->visited_graph[i_next_graph];
					if(num_nn_graph_visited == 0 || i_graph_chain_visited_shortest_loop < i)
					{
						i_nn_graph_visited_shortest_loop = i_nn;
						i_graph_chain_visited_shortest_loop = i;
					}
					num_nn_graph_visited ++;
				}
				if(wrk->is_defect_free[i_next_graph])
				{
					wrk->index_nn_graph_defect_free[num_nn_graph_defect_free] = i_nn;
					num_nn_graph_defect_free ++;
				}
			}
	
			//stop search if acceptable loop or chain found
			if(num_nn_graph_visited > 0)
			{
				i_nn = i_nn_graph_visited_shortest_loop;
				i_flip_start_graph_chain = i_graph_chain_visited_shortest_loop;
				i_flip_end_graph_chain = num_graph_visited - 1;
				move_to_next_graph(&i_graph, &i_nn, graph_n->index_nn_graph, &num_graph_visited, wrk->visited_graph, wrk->graph_chain,
					&num_site_chain, wrk->site_chain, &i_site_last_passed);

				i_flip_start_site_chain = num_site_chain - (num_graph_visited - i_flip_start_graph_chain) + 1;
				i_flip_end_site_chain = num_site_chain - 1;

				exit_mode = VISITED_FOUND;
				break;
			}
			else if(num_nn_graph_defect_free == 0)
			{
				//stop search if no more defect-free graph
				exit_mode = FAILED;
				break;
			}

			//move to another defect-free graph and continue search
			i_nn = wrk->index_nn_graph_defect_free[dsfmt_genrand_uint32(dsfmt)%num_nn_graph_defect_free];
			move_to_next_graph(&i_graph, &i_nn, graph_n->index_nn_graph, &num_graph_visited, wrk->visited_graph, wrk->graph_chain,
				&num_site_chain, wrk->site_chain, &i_site_last_passed);
		}

		//if acceptable loop found, try update
		*num_site_visited += num_site_chain;
		if(exit_mode == VISITED_FOUND)
		{
			int len_loop=i_flip_end_site_chain-i_flip_start_site_chain+1;
			itmp = loop_flip(dsfmt, loop_flip_method, wrk->easy_axis_site,num_update_fix_loop, beta, *p_replica,
				&(wrk->site_chain[i_flip_start_site_chain]),len_loop,
				ham, idx_cnct_bond, idx_cnctd_site, coord_num, rwork);

			(*num_loop_found) ++;
			if(itmp > 0)
				(*num_loop_flipped) ++;
			(*num_site_flipped) += itmp;
		}

		//clean up work array
		for(i = 0; i <= num_graph_visited - 1; i ++)
		{
			wrk->visited_graph[wrk->graph_chain[i]] = -1;
		}
		i_loop ++;
	}
	//cout << "num_site_visited= " << *num_site_visited << endl;
	*num_try = i_loop;
}
