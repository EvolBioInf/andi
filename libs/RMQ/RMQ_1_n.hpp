#ifndef _RMQ_1_n_hpp_
#define _RMQ_1_n_hpp_

#include "RMQ.hpp"
#include <stdlib.h>
#include <iostream>
using namespace std;


/* Implements the most naive <O(1),O(n)>-method to compute RMQ.
 */
class RMQ_1_n : public RMQ {
public:
	// returns RMQ(i,j)
	virtual DTidx query(DTidx, DTidx);

	RMQ_1_n(DT* a, DTidx n);

	~RMQ_1_n();

protected:
	// array
	DT* a;

	// size
	DTidx n;
};

#endif
