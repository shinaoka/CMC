#ifndef PHYS_VALUE_RECORDER_H
#define PHYS_VALUE_RECORDER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <valarray>
#include "phys_value.h"
#include "ham.h"
#include "replica.h"
#include "util.h"

using namespace std;

#define P_NUM_BUILT_IN_PHYS_VALUE 12

class PhysValueCalc
{
public:
	PhysValueCalc()
	{
		ptr_replica = NULL;
		ptr_ham = NULL;
	}

	PhysValueCalc(int num_replica, Ham* ptr_ham)
	{
		set_num_replica(num_replica);
		ptr_replica = new Replica*[num_replica];
		this->ptr_ham = ptr_ham;
	}

	~PhysValueCalc()
	{
		delete[] ptr_replica;
	}

	void add_ptr_replica(int i_replica, Replica *replica)
	{
		if(i_replica < 0 || i_replica >= num_replica)
			throw "PhysValueCalc::add_ptr_replica(): out of range (i_replica)\n";

		ptr_replica[i_replica] = replica;
	}

	inline Ham* get_ptr_ham()
	{
		return ptr_ham;
	}

	inline Replica* get_ptr_replica(int i_replica)
	{
		return ptr_replica[i_replica];
	}

	inline void set_num_replica(int num_replica)
	{
		this->num_replica = num_replica;
	}

	inline int get_num_replica()
	{
		return num_replica;
	}

	virtual int get_num_data_set() = 0;
	virtual int get_num_phys_value() = 0;
	virtual void calc(vector<double>& phys_value, vector<double>& energy) = 0;

private:
	int num_replica;
	Replica **ptr_replica;
	Ham *ptr_ham;
};
typedef class PhysValueCalc PhysValueCalc;

class BuiltInPhysValueCalc : public PhysValueCalc
{
public:
	BuiltInPhysValueCalc(int num_replica, Ham* ptr_ham, int num_sbl_mag, vector<complex<double > > &sbl_mag_coeff, vector<CROSS_PRODUCT_TYPE>& cross_product, bool abs_phys) : PhysValueCalc(num_replica, ptr_ham) {
		int Ns = ptr_ham->Ns;

		if(sbl_mag_coeff.size() < num_sbl_mag*ptr_ham->Ns) {
			cout << "Error :: sbl_mag_coeff.size < num_sbl_mag*Ns!" << endl;
			exit(1);
		}
		this->num_sbl_mag = num_sbl_mag;
		this->sbl_mag_coeff.resize(num_sbl_mag*Ns);
		for(int i=0; i<Ns*num_sbl_mag; ++i){
			this->sbl_mag_coeff[i] = sbl_mag_coeff[i];
		}
		num_phys_value = P_NUM_BUILT_IN_PHYS_VALUE+4*num_sbl_mag+4*cross_product.size();
		num_built_in_phys_value = P_NUM_BUILT_IN_PHYS_VALUE;
		this->abs_phys = abs_phys;

                this->cross_product = cross_product;
	}

	virtual int get_num_data_set()
	{
		return get_num_replica();
	}

	virtual int get_num_phys_value()
	{
		return num_phys_value;
	}

	virtual int get_num_sbl_mag()
	{
		return num_sbl_mag;
	}

	virtual int get_num_built_in_phys_value()
	{
		return P_NUM_BUILT_IN_PHYS_VALUE;
	}

	virtual void calc(vector<double>& phys_value, vector<double>& energy);

private:
	int num_phys_value, num_built_in_phys_value, num_sbl_mag;
	bool abs_phys;
	vector<complex<double> > sbl_mag_coeff;
        vector<CROSS_PRODUCT_TYPE> cross_product;
};

class SGCalc : public PhysValueCalc
{
public:
	SGCalc(int num_replica, int num_seperated_replica_set, Ham* ptr_ham) : PhysValueCalc(num_replica, ptr_ham)
	{
		num_data_set = num_replica*(num_seperated_replica_set - 1)/2;
		this->num_seperated_replica_set = num_seperated_replica_set;
		this->num_replica_per_set = num_replica/num_seperated_replica_set;
	}

	virtual int get_num_data_set()
	{
		return num_data_set;
	}

	virtual int get_num_phys_value()
	{
		return 1;
	}

	virtual void calc(vector<double>& phys_value, vector<double>& energy);
private:
	int num_data_set, num_seperated_replica_set, num_replica_per_set;
};
//typedef class SGPhysValueCalc SGPhysValueCalc;

class SSCCalc : public PhysValueCalc
{
public:
	SSCCalc(int num_replica, Ham* ptr_ham, int num_ss_pair, int **ss_pair) : PhysValueCalc(num_replica, ptr_ham)
	{
		this->num_ss_pair = num_ss_pair;
		this->ss_pair = ss_pair;
	}

	virtual int get_num_data_set()
	{
		return get_num_replica();
	}

	virtual int get_num_phys_value()
	{
		return num_ss_pair;
	}

	virtual void calc(vector<double>& phys_value, vector<double>& energy);
private:
	int num_ss_pair, **ss_pair;
};

class PhysValueRecorder
{
public:
	PhysValueRecorder(PhysValueCalc *phys_value_calc, int num_power, int max_num_sample, int num_sample_disregarded, int reweightable)
	{
		this->phys_value_calc = phys_value_calc;
		this->num_power = num_power;
		this->max_num_sample = max_num_sample;
		this->num_sample_disregarded = num_sample_disregarded;
		this->reweightable = reweightable;
		num_sample = 0;

		num_phys_value = this->phys_value_calc->get_num_phys_value();
		num_data_set = this->phys_value_calc->get_num_data_set();

		if(reweightable){
			phys_value.resize(max_num_sample*num_data_set*num_phys_value);
			energy.resize(max_num_sample*num_data_set);
		}

		buffer.resize(num_data_set*num_phys_value);
		buffer_e.resize(num_data_set);

		sum_phys_value.resize(num_power);
		for(int i_power=0; i_power < num_power; ++i_power)
			sum_phys_value[i_power].resize(num_phys_value);
	}

	PhysValueRecorder()
	{
		this->num_data_set = 0;
		num_sample = 0;
		phys_value_calc = NULL;
	}

	~PhysValueRecorder()
	{
	}

	PhysValueRecorder(vector<double>::iterator& it_begin, vector<double>::iterator& it_end);

	inline double& energy_at(int i_data_set, int i_sample)
	{
		return energy[i_data_set + i_sample * num_data_set];
	}

	inline double& operator()(int i_data_set, int i_phys_value, int i_sample)
	{
		if(i_data_set < 0 || i_data_set >= num_data_set)
			throw "PhysValueRecorder::operator(): out of range (i_data_set)\n";
		if(i_phys_value < 0 || i_phys_value >= num_phys_value)
			throw "PhysValueRecorder::operator(): out of range (i_phys_value)\n";
		if(i_sample < 0 || i_sample >= max_num_sample)
		{
			cout << "i_sample=" << i_sample << " max_num_sample=" << max_num_sample << endl;
			throw "PhysValueRecorder::operator(): out of range (i_sample)\n";
		}

		if(phys_value.size() <= i_phys_value + i_data_set * num_phys_value + i_sample * num_data_set * num_phys_value)
		{
		 	cout << phys_value.size() << " " << i_phys_value + i_data_set * num_phys_value + i_sample * num_data_set * num_phys_value << endl;
			throw "PhysValueRecorder::operator(): out of range\n";
		}
		return phys_value[i_phys_value + i_data_set * num_phys_value + i_sample * num_data_set * num_phys_value];
	}

	int get_num_sample() {return num_sample;}
	int get_num_data_set() {return num_data_set;}
	int get_num_phys_value() {return num_phys_value;}
	int get_num_power() {return num_power;}

	void sampling();
	//void average(int num_sample_disregarded, vector<vector<double> >& average);
	void average(vector<vector<double> >& average);
	void reweight(int num_sample_disregarded, vector<vector<double> >& average, double dE);
	void serialize(vector<double>& double_array);
	void dump(ofstream& stream);
	static PhysValueRecorder& load(ifstream& stream);
private:
	int num_power, num_data_set, num_phys_value;
	int num_sample, max_num_sample, reweightable, num_sample_disregarded;
	vector<double> phys_value, energy;
	vector<vector<double> > sum_phys_value;
	PhysValueCalc *phys_value_calc;
	vector<double> buffer, buffer_e;
};
typedef class PhysValueRecorder PhysValueRecorder;

#endif
