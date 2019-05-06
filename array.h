#include <boost/multi_array.hpp>

#ifndef ARRAY_H
#define ARRAY_H

typedef boost::multi_array<int, 1> int1;
typedef boost::multi_array<int, 2> int2;
typedef boost::multi_array<int, 3> int3;

typedef boost::multi_array<double, 1> double1;
typedef boost::multi_array<double, 2> double2;
typedef boost::multi_array<double, 3> double3;
#endif
