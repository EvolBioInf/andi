#ifndef _RMQ_succinct_hpp_
#define _RMQ_succinct_hpp_

//#define MEM_COUNT

#include "RMQ.hpp"
#include <stdlib.h>
#include <limits.h>
#include <iostream>
using namespace std;

typedef unsigned char DTsucc;
typedef unsigned short DTsucc2;

class RMQ_succinct : public RMQ {
public:
	// liefert RMQ[i,j]
	virtual DTidx query(DTidx, DTidx);

	RMQ_succinct(DT* a, DTidx n);

	~RMQ_succinct();

protected:
	// array
	DT *a;

	// size of array a
	DTidx n;

	// table M for the out-of-block queries (contains indices of block-minima)
	DTsucc** M;

	// because M just stores offsets (rel. to start of block), this method
	// re-calculates the true index:
	inline DTidx m(DTidx k, DTidx block) { return M[k][block]+(block*sprime); }

	// depth of table M:
	DTidx M_depth;

	// table M' for superblock-queries (contains indices of block-minima)
	DTidx** Mprime;

	// depth of table M':
	DTidx Mprime_depth;

	// type of blocks
	DTsucc2 *type;

	// precomputed in-block queries
	DTsucc** Prec;

	// microblock size
	DTidx s;

	// block size
	DTidx sprime;

	// superblock size
	DTidx sprimeprime;

	// number of blocks (always n/sprime)
	DTidx nb;

	// number of superblocks (always n/sprimeprime)
	DTidx nsb;

	// number of microblocks (always n/s)
	DTidx nmb;

	// return microblock-number of entry i:
	inline DTidx microblock(DTidx i) { return i/s; }

	// return block-number of entry i:
	inline DTidx block(DTidx i) { return i/sprime; }

	// return superblock-number of entry i:
	inline DTidx superblock(DTidx i) { return i/sprimeprime; }

	// precomputed Catalan triangle (17 is enough for 64bit computing):
	static const DTidx Catalan[17][17];

	// minus infinity (change for 64bit version)
	static const DT minus_infinity;

 	// stuff for clearing the least significant x bits (change for 64-bit computing)
	static const DTsucc HighestBitsSet[8];
	virtual DTsucc clearbits(DTsucc, DTidx);

	// Least Significant Bits for 8-bit-numbers:
	static const char LSBTable256[256];

    // return least signigicant bit in constant time (change for 64bit version)
	virtual DTidx lsb(DTsucc);

	// the following stuff is for fast base 2 logarithms:
	// (currently only implemented for 32 bit numbers)
	static const char LogTable256[256];

	virtual DTidx log2fast(DTidx);
};

#endif
