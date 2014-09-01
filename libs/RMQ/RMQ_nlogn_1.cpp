#include "RMQ_nlogn_1.hpp"

const char RMQ_nlogn_1::LogTable256[256] = 
	{
		0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
	};

DTidx RMQ_nlogn_1::log2fast(DTidx v) {
	DTidx c = 0;          // c will be lg(v)
	register DTidx t, tt; // temporaries

	if (tt = v >> 16)
		c = (t = v >> 24) ? 24 + LogTable256[t] : 16 + LogTable256[tt & 0xFF];
	else 
		c = (t = v >> 8) ? 8 + LogTable256[t] : LogTable256[v];
	return c;
}

DTidx RMQ_nlogn_1::query(DTidx i, DTidx j) {
	if (j-i == 0) return j;
	if (j-i == 1) return M[0][i];
	DTidx k = log2fast(j-i);
	DTidx twotothek = 1 << k; // 2^k
	return a[c[M[k-1][i]]] <= a[c[M[k-1][j+1-twotothek]]] ? M[k-1][i] : M[k-1][j+1-twotothek];
}

/**
 * Standard Constructor. a[c[0]],...,a[c[n-1]] is the array to be prepared for RMQ.
 * n is the size of the index array c.
 */
RMQ_nlogn_1::RMQ_nlogn_1(DT* a, DTidx* c, DTidx n) {
	this->a = a;
	this->c = c;
	this->n = n;
	depth = log2fast(n); // table depth

	// allocate space for table:
	M = new DTidx*[depth];
	for (DTidx i = 0; i < depth; i++)
		M[i] = new DTidx[n];

	// fill table:
	for (DTidx i = 0; i < n-1; i++) // fill first row
		M[0][i] = a[c[i]] <= a[c[i+1]] ? i : i+1;
	if (depth > 0) M[0][n-1] = n-1;          // fill overhang in first row

	DTidx dist = 1; // always 2^j
	for (DTidx j = 1; j < depth; j++) {
		dist *= 2;
		for (DTidx i = 0; i < n - dist; i++) // fill jth row
			M[j][i] = a[c[M[j-1][i]]] <= a[c[M[j-1][i+dist]]] ? M[j-1][i] : M[j-1][i+dist];
		for (DTidx i = n - dist; i < n; i++) M[j][i] = M[j-1][i]; // overhang
	}
}

/**
 * Destructor. Deletes allocated space.
 */
RMQ_nlogn_1::~RMQ_nlogn_1() {
	for (DTidx i = 0; i < depth; i++)
		delete[] M[i];
	delete[] M;
}
