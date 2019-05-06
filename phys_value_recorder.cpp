#include "phys_value_recorder.h"

void PhysValueRecorder::sampling()
{
	double rtmp;

	if(num_sample == max_num_sample)
		return;

	phys_value_calc->calc(buffer, buffer_e);
	for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set)
	{
		if(num_sample+1 > num_sample_disregarded){
			for(int i_phys_value = 0; i_phys_value < num_phys_value; ++ i_phys_value)
			{
				double rtmp = buffer[i_phys_value + num_phys_value * i_data_set];
				double rtmp2 = rtmp;
				for(int i_power=0; i_power < num_power; ++i_power){
					sum_phys_value[i_power][i_phys_value] += rtmp;
					rtmp *= rtmp2;
				}
			}
		}

		if(reweightable){
			(*this).energy_at(i_data_set, num_sample) = buffer_e[i_data_set];
			for(int i_phys_value = 0; i_phys_value < num_phys_value; ++ i_phys_value){
				(*this)(i_data_set, i_phys_value, num_sample) = buffer[i_phys_value + num_phys_value * i_data_set];
			}
		}
	}
	++ num_sample;
}

void PhysValueRecorder::average(vector<vector<double> >& average)
{
	average.resize(num_power);
	for(int i_power = 0; i_power < num_power ; ++ i_power)
		average[i_power].resize(num_phys_value);

	for(int i_phys_value = 0; i_phys_value < num_phys_value; ++ i_phys_value)
	{
		for(int i_power=0; i_power < num_power; ++i_power){
			average[i_power][i_phys_value] = sum_phys_value[i_power][i_phys_value]*
				(1.0/num_data_set)*(1.0/(num_sample-num_sample_disregarded));
		}
	}
}

void PhysValueRecorder::reweight(int num_sample_disregarded, vector<vector<double> >& average, double dbeta)
{
	vector<vector<double> > sum_phys_value;
	double rtmp, rtmp2;

	if(num_sample <= num_sample_disregarded)
	{
		cout << "num_sample=" << num_sample << endl;
		throw "Too much num_sample_disregarded!\n";
	}
	if(num_sample == 0)
		throw "There is no available sample!\n";

	double ave_E = 0.0;
	for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set){
		for(int i_sample = num_sample_disregarded; i_sample < num_sample; ++ i_sample)
			ave_E += (*this).energy_at(i_data_set, i_sample);
	}
	ave_E /= 1.0*num_data_set*(num_sample - num_sample_disregarded);

	rtmp2 = 0.0;
	for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set)
	{
		for(int i_sample = num_sample_disregarded; i_sample < num_sample; ++ i_sample)
		{
			//rtmp2 += exp(-dbeta*(*this).energy_at(i_data_set, i_sample));
			rtmp2 += exp(-dbeta*((*this).energy_at(i_data_set, i_sample)-ave_E));
		}
	}
	rtmp2 /= 1.0*num_data_set*(num_sample - num_sample_disregarded);

	sum_phys_value.resize(num_power);
	for(int i_power = 0; i_power < num_power ; ++ i_power)
		sum_phys_value[i_power].resize(num_phys_value, 0.0);
	
	for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set){
		for(int i_phys_value = 0; i_phys_value < num_phys_value ; ++ i_phys_value){
			for(int i_sample = num_sample_disregarded; i_sample < num_sample; ++ i_sample){
				double t_phys = (*this)(i_data_set, i_phys_value, i_sample);
				double rtmp = t_phys;
				//double t_exp = exp(-dbeta*(*this).energy_at(i_data_set, i_sample));
				double t_exp = exp(-dbeta*((*this).energy_at(i_data_set, i_sample)-ave_E));
				for(int i_power = 0; i_power < num_power ; ++ i_power){
					sum_phys_value[i_power][i_phys_value] += t_exp * t_phys;
					t_phys *= rtmp;
				}
			}
		}
	}
	for(int i_power = 0; i_power < num_power; ++ i_power){
		for(int i_phys_value = 0; i_phys_value < num_phys_value ; ++ i_phys_value){
			sum_phys_value[i_power][i_phys_value] *= 1.0/((num_sample - num_sample_disregarded)*rtmp2);
		}
	}

	average.resize(num_power);
	for(int i_power = 0; i_power < num_power ; ++ i_power){
		average[i_power].resize(num_phys_value);
		for(int i_phys_value = 0; i_phys_value < num_phys_value; ++ i_phys_value){
			average[i_power][i_phys_value] = sum_phys_value[i_power][i_phys_value]/num_data_set;
		}
	}
}

void PhysValueRecorder::serialize(vector<double>& double_array)
{
	int t_num_sample;

	double_array.resize(0);
	double_array.push_back((double) num_power);
	double_array.push_back((double) num_data_set);
	double_array.push_back((double) num_phys_value);
	double_array.push_back((double) num_sample);
	double_array.push_back((double) max_num_sample);
	double_array.push_back((double) num_sample_disregarded);
	double_array.push_back((double) reweightable);
	
	if(reweightable){
		for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set){
			for(int i_sample = 0; i_sample < num_sample; ++ i_sample)
				double_array.push_back((*this).energy_at(i_data_set, i_sample));
		}
	
		for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set){
			for(int i_phys_value = 0; i_phys_value < num_phys_value; ++ i_phys_value){
				for(int i_sample = 0; i_sample < num_sample; ++ i_sample)
					double_array.push_back((*this)(i_data_set, i_phys_value, i_sample));
			}
		}
	}

	for(int i_power=0; i_power < num_power; ++i_power){
		for(int i_phys_value = 0; i_phys_value < num_phys_value; ++i_phys_value)
			double_array.push_back(sum_phys_value[i_power][i_phys_value]);
	}
}

PhysValueRecorder::PhysValueRecorder(vector<double>::iterator& it_begin,
	vector<double>::iterator& it_end)
{
	vector<double>::iterator it = it_begin;

	if(distance(it_begin, it_end) < 1)
		throw "PhysValueRecorder:: too short input vector (header)\n";

	num_power = (int) *(it++);
	num_data_set = (int) *(it++);
	num_phys_value = (int) *(it++);
	num_sample = (int) *(it++);
	max_num_sample = (int) *(it++);
	num_sample_disregarded = (int) *(it++);
	reweightable = (int) *(it++);

	if(reweightable){
		energy.resize(num_data_set*max_num_sample);
		for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set){
			for(int i_sample = 0; i_sample < num_sample; ++ i_sample)
				(*this).energy_at(i_data_set, i_sample) = *(it++);
		}
	
		phys_value.resize(num_data_set*num_phys_value*max_num_sample);
		for(int i_data_set = 0; i_data_set < num_data_set; ++ i_data_set){
			for(int i_phys_value = 0; i_phys_value < num_phys_value; ++ i_phys_value){
				for(int i_sample = 0; i_sample < max_num_sample; ++ i_sample)
					(*this)(i_data_set, i_phys_value, i_sample) = *(it++);
			}
		}
	}

	sum_phys_value.resize(num_power);
	for(int i_power=0; i_power < num_power; ++i_power){
		sum_phys_value[i_power].resize(num_phys_value);

		for(int i_phys_value = 0; i_phys_value < num_phys_value; ++i_phys_value)
			sum_phys_value[i_power][i_phys_value] = *(it++); 
	}
}

PhysValueRecorder& PhysValueRecorder::load(ifstream& stream)
{
	vector<double> buffer;
	string str;
	int num_data;

	stream.unsetf(ios::unitbuf);

	stream >> num_data;
	buffer.resize(num_data);
	for(int i=0; i < num_data; i++)
		stream >> buffer[i];
	vector<double>::iterator it_begin = buffer.begin();
	vector<double>::iterator it_end = buffer.end();
	return *(new PhysValueRecorder(it_begin, it_end));
}

void PhysValueRecorder::dump(ofstream& stream)
{
	vector<double> buffer;

	serialize(buffer);
	stream.setf(ios::showpos);
	stream.setf(ios::scientific);
	stream.unsetf(ios::unitbuf);
	vector<double>::iterator it = buffer.begin();
	stream << buffer.size() << '\n';
	while(it != buffer.end())
	{
		stream << *it << '\n';
		++ it;
	}
}

void BuiltInPhysValueCalc::calc(vector<double>& buffer, vector<double>& energy)
{
	//double t_phys_value[NUM_BUILT_IN_PHYS_VALUE-1];
	vector<double> t_phys_value;
	vector<complex<double > > t_sbl_mag;
	vector<complex<double > > t_cross_product;

	int num_replica;
	Ham *ptr_ham;
	Replica *ptr_replica;

	num_replica = get_num_replica();

	if(buffer.size() < num_replica*num_phys_value)
		buffer.resize(num_replica*num_phys_value);
	if(energy.size() < num_replica)
		energy.resize(num_replica);
	ptr_ham = get_ptr_ham();
	t_phys_value.resize(num_built_in_phys_value);
	t_sbl_mag.resize(3*num_sbl_mag);
	t_cross_product.resize(3*cross_product.size());

	for(int i_replica = 0; i_replica < num_replica; ++ i_replica)
	{
		ptr_replica = get_ptr_replica(i_replica);
		buffer[0 + i_replica*num_phys_value] = calc_E(ptr_ham, ptr_replica->config);
		energy[i_replica] = buffer[0 + i_replica*num_phys_value];
		calc_phys_value(ptr_ham->Ns, ptr_replica->config, &t_phys_value[0], ptr_ham->num_two_body_int, ptr_ham->two_body_int_pair,
			num_sbl_mag, sbl_mag_coeff, t_sbl_mag,
			cross_product, t_cross_product,
                        this->abs_phys);
		for(int i=1; i<num_built_in_phys_value; ++i)
			buffer[i + i_replica*num_phys_value] = t_phys_value[i-1];

		for(int i_sbl_mag=0; i_sbl_mag<num_sbl_mag; ++i_sbl_mag){
			double rsum=0.0;
			for(int ix=0; ix<3; ++ix){
				double rtmp = (t_sbl_mag[ix+3*i_sbl_mag].real());
				if (abs_phys) rtmp = fabs(rtmp);
				buffer[ix+4*i_sbl_mag+num_built_in_phys_value+i_replica*num_phys_value] = rtmp;
                                //to be consistent with Andrea's definition
                                //if (ix!=2) rsum += rtmp*rtmp;
                                rsum += rtmp*rtmp;
			}
                        //to be consistent with Andrea's definition
			//buffer[3+4*i_sbl_mag+num_built_in_phys_value+i_replica*num_phys_value] = rsum;
			buffer[3+4*i_sbl_mag+num_built_in_phys_value+i_replica*num_phys_value] = sqrt(rsum);
		}
		for(int icr=0; icr<cross_product.size(); ++icr){
			double rsum=0.0;
			for(int ix=0; ix<3; ++ix){
				double rtmp = (t_cross_product[ix+3*icr].real());
				if (abs_phys) rtmp = fabs(rtmp);
				buffer.at(ix+4*icr+4*num_sbl_mag+num_built_in_phys_value+i_replica*num_phys_value) = rtmp;
				rsum += rtmp*rtmp;
			}
			buffer.at(3+4*icr+4*num_sbl_mag+num_built_in_phys_value+i_replica*num_phys_value) = rsum;
		}
	}
}

void SGCalc::calc(vector<double>& buffer, vector<double>& energy)
{
	int i_data_set = 0, num_data_set;
	Ham *ptr_ham = get_ptr_ham();
	Replica *ptr_replica1, *ptr_replica2;
	double q2, q4;

	num_data_set = get_num_data_set();
	if(buffer.size() < num_data_set)
		buffer.resize(num_data_set);
	if(energy.size() < num_data_set)
		energy.resize(num_data_set);

	for (int i_replica_set = 0; i_replica_set <= num_seperated_replica_set - 1; i_replica_set++)
	{
		for (int i_replica_set2 = i_replica_set + 1; i_replica_set2 <= num_seperated_replica_set - 1; i_replica_set2++)
		{
			for (int i_replica_reduced = 0; i_replica_reduced <= num_replica_per_set - 1; i_replica_reduced++)
			{
				ptr_replica1 = get_ptr_replica(i_replica_reduced + num_replica_per_set*i_replica_set);
				ptr_replica2 = get_ptr_replica(i_replica_reduced + num_replica_per_set*i_replica_set2);
				calc_squared_q
				(
					ptr_ham->Ns, ptr_replica1, ptr_replica2, &q2, &q4
				);
				buffer[i_data_set] = q2;
				energy[i_data_set] = ptr_replica1->E + ptr_replica2->E;
				++ i_data_set;
			}
		}
	}
}

void SSCCalc::calc(vector<double>& buffer, vector<double>& energy)
{
	int num_replica = get_num_replica();
	if(buffer.size() < num_replica*num_ss_pair)
		buffer.resize(num_replica*num_ss_pair);
	if(energy.size() < num_replica)
		energy.resize(num_replica);

	for(int i_replica = 0; i_replica < num_replica; ++ i_replica){
		Replica *ptr_replica = get_ptr_replica(i_replica);
		double **config = ptr_replica->config;

		energy[i_replica] = ptr_replica->E;
		for(int i = 0; i < num_ss_pair; ++ i){
			int i_site0 = ss_pair[i][0];
			int i_site1 = ss_pair[i][1];
			buffer[i + i_replica*num_ss_pair] = dot_product(config[i_site0], config[i_site1]);
		}
	}
}
