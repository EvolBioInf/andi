#include "RMQ_n_1_improved.hpp"

const DTidx RMQ_n_1_improved::Catalan[17][17] = {
	{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16},
	{0,0,2,5,9,14,20,27,35,44,54,65,77,90,104,119,135},
	{0,0,0,5,14,28,48,75,110,154,208,273,350,440,544,663,798},
	{0,0,0,0,14,42,90,165,275,429,637,910,1260,1700,2244,2907,3705},
	{0,0,0,0,0,42,132,297,572,1001,1638,2548,3808,5508,7752,10659,14364},
	{0,0,0,0,0,0,132,429,1001,2002,3640,6188,9996,15504,23256,33915,48279},
	{0,0,0,0,0,0,0,429,1430,3432,7072,13260,23256,38760,62016,95931,144210},
	{0,0,0,0,0,0,0,0,1430,4862,11934,25194,48450,87210,149226,245157,389367},
	{0,0,0,0,0,0,0,0,0,4862,16796,41990,90440,177650,326876,572033,961400},
	{0,0,0,0,0,0,0,0,0,0,16796,58786,149226,326876,653752,1225785,2187185},
	{0,0,0,0,0,0,0,0,0,0,0,58786,208012,534888,1188640,2414425,4601610},
	{0,0,0,0,0,0,0,0,0,0,0,0,208012,742900,1931540,4345965,8947575},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,742900,2674440,7020405,15967980},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,2674440,9694845,25662825},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9694845,35357670},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,35357670}
};

const DT RMQ_n_1_improved::minus_infinity = -2147483648;

DTidx RMQ_n_1_improved::query(DTidx i, DTidx j) {
	DTidx block_i = block(i);       // block where i is
	DTidx block_j = block(j);       // j's block
	DTidx i_pos = i - (block_i*bs); // position of i in block_i
	DTidx j_pos = j - (block_j*bs); // pos. of j in block_j
	DTidx min, min_i, min_j; 		// min: to be returned

	if (block_j == block_i)
		min = (block_i*bs) + Prec[type[block_i]][i_pos*(bs-1)+j_pos-((i_pos-1)*i_pos/2)];
	else {
		min_i = (block_i*bs) + Prec[type[block_i]][(i_pos+1)*(bs-1)-((i_pos-1)*i_pos/2)];
		min_j = (block_j*bs) + Prec[type[block_j]][j_pos];
		min = a[min_i] < a[min_j] ? min_i : min_j;
		if (block_j > block_i + 1) { // out-of-block query:
			min_i = c[RMQ_ST->query(block_i+1, block_j-1)];
			if (a[min_i] < a[min]) min = min_i;
		}
	}

	return min;
}

/**
 * Standard Constructor. a is the array to be prepared for RMQ.
 * n is the size of the array.
 */
RMQ_n_1_improved::RMQ_n_1_improved(DT* a, DTidx n) {
	this->a = a;
	this->n = n;
	bs = (DTidx) ceil(log2(n) / 4); // calculate block size
	nb = block(n-1)+1;              // number of blocks
	qpb = bs * (bs+1) / 2;          // # of different queries per block

	// allocate and fill arrays for out-of-block-queries;
	// combined with type-calculation for the blocks:
	c = new DTidx[nb];
	type = new DTidx[nb];
	Prec = new DTidx*[Catalan[bs][bs]];
	for (DTidx i = 0; i < Catalan[bs][bs]; i++) {
		Prec[i] = new DTidx[qpb];
		Prec[i][0] = bs; // init with impossible value
	}
	DT* rp = new DT[bs+1];  // rp: rightmost path in Cart. tree

	DTidx min;   // position of minimum in current block
	DTidx z=0;   // index in array a
	DTidx end;   // end of current block
	DTidx start; // start of current block
	DTidx q;     // position in Catalan triangle
	DTidx p;     // --------- " ----------------
	rp[0] = minus_infinity; // stopper (minus infinity)

	for (DTidx i = 0; i < nb; i++) { // step through blocks
		start = z;            // init start
		min = start;          // init minimum
		end = start + bs;     // end of block
		if (end > n) end = n; // last block could be smaller than bs!
		q = bs;       // init q
		p = bs-1;     // init p
		type[i] = 0;  // init type (will be increased!)
		rp[1] = a[z]; // init rightmost path

		while (++z < end) {   // step through current block:
			if (a[z] < a[min]) min = z;   // update minimum
			p--;
			while (rp[q-p-1] > a[z]) {
				type[i] += Catalan[p][q]; // update type
				q--;
			}
			rp[q-p] = a[z]; // add last element to rightmost path
		}
		c[i] = min; // store index of minimum

		// precompute in-block-queries for this block (if necessary):
		if (Prec[type[i]][0] == bs) {
			p = 0; // now index in Prec[type[i]]
			for (DTidx j = start; j < end; j++) {
				min = j;
				for (DTidx k = j; k < end; k++) {
					if (a[k] < a[min]) min = k;
					Prec[type[i]][p++] = min-start;
				}
			}
		}
	}
	delete[] rp;

	RMQ_ST = new RMQ_nlogn_1(a, c, nb); // construct sparse table on block-minima
}

/**
 * Destructor. Deletes allocated space.
 */
RMQ_n_1_improved::~RMQ_n_1_improved() {
	delete RMQ_ST;
	delete[] type;
	delete[] c;
	for (DTidx i = 0; i < Catalan[bs][bs]; i++)
		delete[] Prec[i];
	delete[] Prec;
}
