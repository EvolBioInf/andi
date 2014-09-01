#ifndef _RMQ_nlogn_1_hpp_
#define _RMQ_nlogn_1_hpp_

#include "RMQ.hpp"
#include <stdlib.h>
#include <iostream>
using namespace std;


/* Implements the <O(n log n), O(1)>-method for RMQ as described in
 * Bender and Farach's paper, section 3.
 */
class RMQ_nlogn_1 : public RMQ {
public:
	// liefert RMQ[i,j]
	virtual DTidx query(DTidx, DTidx);

	RMQ_nlogn_1(DT* a, DTidx* c, DTidx n);

	~RMQ_nlogn_1();

	// the following stuff is for fast base 2 logarithms:
	// (currently only implemented for 32 bit numbers)
	static const char LogTable256[256];

	virtual DTidx log2fast(DTidx);

protected:
	// array
	DT* a;

	// index array for a:
	DTidx* c;

	// size of c
	DTidx n;

	// depth of table
	DTidx depth;

	// the precomputed table:
	DTidx** M;

};

#endif
