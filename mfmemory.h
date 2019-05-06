/*-------------------------------------------------------------
 * Memmory allocator for multi-dimensional arrays (compatible with Taraha's version)
 *-------------------------------------------------------------
 * Copyright (C) 2010-2012 Hiroshi Shinaoka. All rights reserved.
 *-------------------------------------------------------------*/

/*=================================================================================================*/
#ifndef _INCLUDE_MFMEMORY
#define _INCLUDE_MFMEMORY

template<class T>
class Matrix1D
{
public:
	Matrix1D(int N1)
	{
		array = new T[N1];
		this->N1 = N1;
	}

	~Matrix1D()
	{
		delete(array);
	}

	inline T& operator()(int i)
	{
		return array[i];
	}

	inline T& operator()(int& i)
	{
		return array[i];
	}

private:
	int N1;
	T* array;
};

template<class T>
class Matrix2D
{
public:
	Matrix2D(int N1, int N2)
	{
		array = new T[N1*N2];
		this->N1 = N1;
		this->N2 = N2;
	}

	~Matrix2D()
	{
		delete(array);
	}

	inline T& operator()(int i, int j)
	{
		return array[i*N2 + j];
	}

	inline T& operator()(int& i, int& j)
	{
		return array[i*N2 + j];
	}

private:
	int N1, N2;
	T* array;
};

template<class T>
class Matrix3D
{
public:
	Matrix3D(int N1, int N2, int N3)
	{
		array = new T[N1*N2*N3];
		this->N1 = N1;
		this->N2 = N2;
		this->N3 = N3;
	}

	~Matrix3D()
	{
		delete(array);
	}

	inline T& operator()(int i, int j, int k)
	{
		return array[i*N2*N3 + j*N3 + k];
	}

	inline T& operator()(int& i, int& j, int& k)
	{
		return array[i*N2*N3 + j*N3 + k];
	}

private:
	int N1, N2, N3;
	T* array;
};


/****
	Wrapper macro
****/
/*double type*/
#define d_malloc1(X, N1) \
	malloc1_seq((X), (N1));

#define d_free1(X, N1) \
	free1_seq((X));

#define d_malloc2(X, N1, N2) \
	malloc2_seq((X), (N1), (N2));

#define d_free2(X, N1, N2) \
	free2_seq((X));

#define d_malloc3(X, N1, N2, N3) \
	malloc3_seq((X), (N1), (N2), (N3));

#define d_free3(X, N1, N2, N3) \
	free3_seq((X));

/*integer type*/
#define i_calloc1(X, N1) \
	X = (int*)calloc((N1), sizeof(int));

#define i_malloc1(X, N1) \
	malloc1_seq((X), (N1));

#define i_free1(X, N1) \
	free1_seq((X));

#define i_malloc2(X, N1, N2) \
	malloc2_seq((X), (N1), (N2));

#define i_free2(X, N1, N2) \
	free2_seq((X));

#define i_malloc3(X, N1, N2, N3) \
	malloc3_seq((X), (N1), (N2), (N3));

#define i_free3(X, N1, N2, N3) \
	free3_seq((X));

/****
	Allocate multi-dimensional array sequantially in the memory space
****/
/******** 1D ********/
template<typename T>
inline void malloc1_seq(T*& ptr, int N1)
{
	if(N1==0){
		ptr = 0;
	}else{
		ptr = new T[N1];
	}
}

template<typename T>
inline void free1_seq(T*& ptr)
{
	if(ptr==0)
		return;
	delete[](ptr);
}

/******** 2D ********/
template<typename T>
inline void malloc2_seq(T**& ptr, int N1, int N2)
{
	int i1, i2;

	if(N1*N2==0){
		ptr = 0;
		return;
	}

	ptr = new T*[N1];

	ptr[0] = new T[N1*N2];
	for(i1 = 0; i1 < N1; i1 ++)
	{
		ptr[i1] = ptr[0] + i1*N2;
	}
}

template<typename T>
inline void free2_seq(T**& ptr)
{
	if(ptr==0)
		return;
	delete[](ptr[0]);
	delete[](ptr);
}

/******** 3D ********/
template<typename T>
inline void malloc3_seq(T***& ptr, int N1, int N2, int N3)
{
	int i1, i2;

	if(N1*N2*N3==0){
		ptr = 0;
		return;
	}

	ptr = new T**[N1];

	ptr[0] = new T*[N1*N2];
	for(i1 = 0; i1 < N1; i1 ++)
	{
		ptr[i1] = ptr[0] + i1*N2;
	}

	ptr[0][0] = new T[N1*N2*N3];
	for(i1 = 0; i1 < N1; i1 ++)
	{
		for(i2 = 0; i2 < N2; i2 ++)
		{
			ptr[i1][i2] = ptr[0][0] + i1*N2*N3 + i2*N3;
		}
	}
}

template<typename T>
inline void free3_seq(T***& ptr)
{
	if(ptr==0)
		return;
	delete[](ptr[0][0]);
	delete[](ptr[0]);
	delete[](ptr);
}

#endif
