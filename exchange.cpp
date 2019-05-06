// int mode=0 (send) or 1 (recv)
#ifdef _mpi_use
int exchange_inter_process(dsfmt_t *dsfmt, int Ns, Replica *p_replica, double my_beta, int i_target, int mode, double *send_buffer, double *recv_buffer)
{
	int i_site, itmp, len_buffer;
	double rtmp;

	MPI_Status status;

	len_buffer = 3*Ns+2;

	if(mode == 0)
	{
		send_buffer[0] = p_replica->E;
		send_buffer[1] = my_beta;
		MPI_Send(send_buffer, 2, MPI_DOUBLE, i_target, 1245, MY_MPI_COMM_WORLD);
		MPI_Recv(recv_buffer, 1, MPI_DOUBLE, i_target, 1246, MY_MPI_COMM_WORLD, &status);
	} else {
		MPI_Recv(recv_buffer, 2, MPI_DOUBLE, i_target, 1245, MY_MPI_COMM_WORLD, &status);
		if(exp((p_replica->E-recv_buffer[0])*(my_beta-recv_buffer[1]))>dsfmt_genrand_close_open(dsfmt))
		{
			*send_buffer = 1.0;
		} else {
			*send_buffer = 0.0;
		}
		MPI_Send(send_buffer, 1, MPI_DOUBLE, i_target, 1246, MY_MPI_COMM_WORLD);
	}

	//exchange E and config
	if((mode==0 && *recv_buffer == 1.0) || (mode==1 && *send_buffer == 1.0))
	{
		itmp = 1;
		send_buffer[0] = p_replica->E;
		for(i_site=0; i_site<=Ns-1; i_site++)
		{
			send_buffer[3*i_site+1] = p_replica->config[i_site][0];
			send_buffer[3*i_site+2] = p_replica->config[i_site][1];
			send_buffer[3*i_site+3] = p_replica->config[i_site][2];
		}
		send_buffer[len_buffer-1] = (double) p_replica->id;

		if(mode == 0)
		{
			MPI_Sendrecv(send_buffer, len_buffer, MPI_DOUBLE, i_target, 1247,
				recv_buffer, len_buffer, MPI_DOUBLE, i_target, 1248,
				MY_MPI_COMM_WORLD, &status);
		} else {
			MPI_Sendrecv(send_buffer, len_buffer, MPI_DOUBLE, i_target, 1248,
				recv_buffer, len_buffer, MPI_DOUBLE, i_target, 1247,
				MY_MPI_COMM_WORLD, &status);
		}
		p_replica->E = recv_buffer[0];
		for(i_site=0; i_site<=Ns-1; i_site++)
		{
			p_replica->config[i_site][0] = recv_buffer[3*i_site+1];
			p_replica->config[i_site][1] = recv_buffer[3*i_site+2];
			p_replica->config[i_site][2] = recv_buffer[3*i_site+3];
		}
		p_replica->id = (int) recv_buffer[len_buffer-1];

	} else {
		itmp = 0;
	}
	return itmp;
}
#endif

void exchange(dsfmt_t *dsfmt, int Ns, int num_beta_per_process, Replica *p_replica,
		int num_thread, int my_rank, int num_beta, double *beta_list, int num_inner_process_ex,
		int *ex_counter, double *send_buffer, double *recv_buffer)
{
	int i_target, i_ex, itmp, len_buffer, mode, i_inner_ex;
	int i_site, i_beta;

	//temporary variables for swap
	double **p_config, rtmp;
	int id;

#ifdef _mpi_use
	MPI_Status status;
#endif

	len_buffer = 3*Ns+2;

	//Inner-process exchange
	//for(i_inner_ex=0; i_inner_ex<=num_inner_process_ex-1; i_inner_ex++)
	//{
		for(i_beta=0; i_beta<=num_beta_per_process-2; i_beta++)
		{
			if(exp((beta_list[i_beta+1+num_beta_per_process*my_rank]-beta_list[i_beta+num_beta_per_process*my_rank])
				*(p_replica[i_beta+1].E-p_replica[i_beta].E)
				)>dsfmt_genrand_close_open(dsfmt))
			{
				ex_counter[i_beta] ++;
	
				p_config = p_replica[i_beta].config;
				rtmp = p_replica[i_beta].E;
				id = p_replica[i_beta].id;
	
				p_replica[i_beta].config = p_replica[i_beta+1].config;
				p_replica[i_beta].E = p_replica[i_beta+1].E;
				p_replica[i_beta].id = p_replica[i_beta+1].id;
	
				p_replica[i_beta+1].config = p_config;
 				p_replica[i_beta+1].E = rtmp;
 				p_replica[i_beta+1].id = id;
			}
		}
	//}

	//Inter-process exchange
	#ifdef _mpi_use
	if(num_thread>1)
	{
		for(i_ex=0; i_ex<=1; i_ex++)
		{
			i_target = my_rank + ((int) pow(-1.0, my_rank%2+i_ex+0.0));	
			if(i_target < 0 || i_target > num_thread-1) continue;
			
			mode = my_rank%2;
			if(i_target == my_rank+1)
			{
				ex_counter[num_beta_per_process-1] += exchange_inter_process(dsfmt, Ns, &p_replica[num_beta_per_process-1]
					, beta_list[(my_rank+1)*num_beta_per_process-1], i_target, mode, send_buffer, recv_buffer);
			} else if(i_target == my_rank-1) {
				exchange_inter_process(dsfmt, Ns, &p_replica[0]
					, beta_list[my_rank*num_beta_per_process], i_target, mode, send_buffer, recv_buffer);
			}
		}
	}
	#endif
}

void opt_beta(int num_beta, int num_beta_per_process, int my_rank, double *beta_list, int min_beta, int max_beta, int *ex_counter)
{
	int *ex_rate = new int[num_beta];
	double *new_beta_list = new double[num_beta];

	#ifdef _mpi_use
	{
		MPI_Status status;
		MPI_Gather(ex_counter, num_beta_per_process, MPI_INT, ex_rate, num_beta_per_process, MPI_INT, 0, MY_MPI_COMM_WORLD);
	}
	#else
	for(int i_beta = 0; i_beta < num_beta; ++ i_beta)
		ex_rate[i_beta] = ex_counter[i_beta];	
	#endif

	if(my_rank == 0)
	{
		double rtmp = 0.0;
		for(int i_beta = min_beta; i_beta < max_beta; ++ i_beta)
		{
			rtmp += (beta_list[i_beta+1] - beta_list[i_beta])*(ex_rate[i_beta]+1);
		}

		//cout << "rtmp=" << rtmp << endl;
		for(int i_beta = 0; i_beta < num_beta; ++ i_beta)
			new_beta_list[i_beta] = beta_list[i_beta];
		if(rtmp > 0.0)
		{
			rtmp = (beta_list[max_beta]-beta_list[min_beta])/rtmp;
			new_beta_list[min_beta] = beta_list[min_beta];
			for(int i_beta = min_beta; i_beta < max_beta; ++ i_beta)
				new_beta_list[i_beta+1] = new_beta_list[i_beta] + (beta_list[i_beta+1] - beta_list[i_beta])*(ex_rate[i_beta]+1)*rtmp;
		}
	}

	#ifdef _mpi_use
	MPI_Bcast(new_beta_list, num_beta, MPI_DOUBLE, 0, MY_MPI_COMM_WORLD);
	#endif
	for(int i_beta = 0; i_beta < num_beta; ++ i_beta)
		beta_list[i_beta] = 0.2*new_beta_list[i_beta] + 0.8*beta_list[i_beta];	

	delete[] ex_rate;
	delete[] new_beta_list;
}
