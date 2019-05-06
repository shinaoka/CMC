#ifndef MPI_UTIL_H
#define MPI_UTIL_H
#ifdef _mpi_use

#include <vector>
#include <mpi.h>
#include "phys_value_recorder.h"

extern MPI_Comm MY_MPI_COMM_WORLD;

using namespace std;

template<typename Iterator, typename T>
void tmpl_mpi_gather_phys_value(MPI_Comm comm, Iterator& input_it, Iterator& output_it, int num_beta, int num_process, int my_rank)
{
	int num_beta_per_process = num_beta/num_process;

	vector<double> snd_buffer, rsv_buffer;

	for(int i_beta = 0; i_beta < num_beta_per_process; ++ i_beta)
	{
		(*(input_it + i_beta)).serialize(snd_buffer);
		int len_buffer = snd_buffer.size();
		if(my_rank == 0)
			rsv_buffer.resize(len_buffer*num_process);

		MPI_Gather(&snd_buffer[0], len_buffer, MPI_DOUBLE, &rsv_buffer[0], len_buffer, MPI_DOUBLE, 0, comm);

		if(my_rank == 0)
		{
			for(int i_process = 0; i_process < num_process; ++ i_process)
			{
				try
				{
					vector<double>::iterator it_begin = rsv_buffer.begin() + len_buffer * i_process;
					vector<double>::iterator it_end = rsv_buffer.begin() + len_buffer * (i_process + 1) - 2;
					*(output_it + i_beta + num_beta_per_process * i_process) = *(new T(it_begin, it_end));
				}
				catch(const char* str)
				{
					cout << str;
					exit(1);
				}
			}
		}
	}
}


#endif
#endif
