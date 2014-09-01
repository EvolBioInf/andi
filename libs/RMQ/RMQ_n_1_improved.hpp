#ifndef _RMQ_n_1_improved_hpp_
#define _RMQ_n_1_improved_hpp_

#include "RMQ.hpp"
#include "RMQ_nlogn_1.hpp"
#include <stdlib.h>
#include <iostream>
using namespace std;

class RMQ_n_1_improved : public RMQ {
public:
	// liefert RMQ[i,j]
	virtual DTidx query(DTidx, DTidx);

	RMQ_n_1_improved(DT* a, DTidx n);

	~RMQ_n_1_improved();

protected:
	// array
	DT *a;

	// index array for the out-of-block queries (contains indices of block-minima)
	DTidx *c;

	// type of blocks
	DTidx *type;

	// precomputed in-block queries
	DTidx** Prec;

	// number of different queries per block = bs*(bs+1)/2
	DTidx qpb;

	// size of array a
	DTidx n;

	// block size
	DTidx bs;

	// number of blocks (always n/bs)
	DTidx nb;

	// return block of entry i:
	inline DTidx block(DTidx i) { return i/bs; }

	// nlogn_1-Algo for out-of-block-block-queries:
	RMQ_nlogn_1* RMQ_ST;

	// precomputed Catalan triangle (17 is enough for 64bit computing):
	static const DTidx Catalan[17][17];

	// minus infinity (change for 64bit version)
	static const DT minus_infinity;
};

#endif
