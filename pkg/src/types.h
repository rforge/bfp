#ifndef TYPES_H_
#define TYPES_H_

#include <set>
#include <vector>
#include <map>

#include <R.h>
#include <newmat.h>

// common type defs:

// the machine precision
static const double EPS = sqrt(DOUBLE_EPS);



// here is a double vector
typedef std::vector<double> DoubleVector;

// and a long double vector
typedef std::vector<long double> LongDoubleVector;




// general signed integer
typedef int Int;

// general unsigned integer
typedef unsigned int PosInt;

// general large unsigned integer
typedef unsigned long long int PosLargeInt;



// here is an integer vector
typedef std::vector<Int> IntVector;

// and an unsigned int vector
typedef std::vector<PosInt> PosIntVector;



// a set of ints
typedef std::set<Int> IntSet;

// a set of unsigned ints
typedef std::set<PosInt> PosIntSet;




// one power set
typedef std::multiset<int> Powers;

// vector of power sets
typedef std::vector<Powers>  PowersVector;




// index type
typedef LongDoubleVector::size_type Index;

// and derived structures
typedef std::set<Index> IndexSet;



// derivatives from matrix types
typedef std::vector<std::vector<ColumnVector> > ColumnVectorArray;







#endif /* TYPES_H_ */
