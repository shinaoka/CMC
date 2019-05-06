#ifndef LOOP_MOVE_H
#define LOOP_MOVE_H

#include <iostream>
#include "dSFMT.h"
#include "ham.h"
#include "replica.h"
#include "mfmemory.h"
#include "phys_value.h"
#include "stdlib.h"
#include "math.h"
#include "util.h"
#include "update.h"

struct GraphNetwork
{
		int num_graph;
		int max_num_nn_graph;
		int *num_nn_graph;
		int ***index_nn_graph;

		int *num_site_graph;
		int max_num_site_graph;
		int **index_site_graph;

		int **bond_graph;

		int **easy_axis_phase; 
};
typedef struct GraphNetwork GraphNetwork;

struct LoopMoveWorkArrays
{
	//work arrays (size == ham->Ns)
	int *visited_graph, *site_chain, *site_mask;
	//work arrays (size == (ham->Ns, 3))
	double **easy_axis_site;
	//work arrays (size == num_graph)
	int *is_defect_free, *index_graph_defect_free, *index_defect_free_graph, *graph_chain;
	//work arrays (size == max_num_nn_graph)
	int *index_nn_graph_visited, *index_nn_graph_correctable_defect, *index_nn_graph_defect_free;
};
typedef struct LoopMoveWorkArrays LoopMoveWorkArrays;

void construct_complete_graph_network(GraphNetwork *graph_n, dsfmt_t *dsfmt, Ham *ham,
	int **idx_cnct_bond, int **idx_cnctd_site, int *num_nn_two_body_int);
void free_complete_graph_network(GraphNetwork *graph_n);

void allocate_loop_move_work_arrays(LoopMoveWorkArrays *wrk, int Ns, int nu_graph, int max_num_nn_graph);
void free_loop_move_work_arrays(LoopMoveWorkArrays *wrk, int Ns, int nu_graph, int max_num_nn_graph);

void loop_move_ISING(dsfmt_t *dsfmt, Replica *p_replica, double beta,
		Ham *ham, int **idx_cnct_bond, int **idx_cnctd_site,
		int *num_nn_two_body_int, double **ia,
		int num_loop, int num_update_fix_loop, int num_graph_ising_axis, int loop_flip_method,
		GraphNetwork *graph_n, LoopMoveWorkArrays *wrk,
		int *num_site_visited, int *num_try, int *num_loop_found, int *num_loop_flipped, int *num_site_flipped);
#endif
