#ifndef _ESA_H_
#define _ESA_H_

#include <divsufsort.h>
#include <RMQ.hpp>


typedef struct {
	const char *S;
	saidx_t *SA;
	saidx_t *ISA;
	saidx_t *LCP;
	saidx_t len;
	RMQ *rmq_lcp;
} esa_t;

typedef struct {
	saidx_t i;
	saidx_t j;
} interval;

typedef struct {
	saidx_t l, i, j;
} lcp_inter_t;


int compute_SA( esa_t *c);
int compute_LCP( esa_t *c);
int compute_LCP_PHI( esa_t *c);
saidx_t longestMatch( const esa_t *C, const char *query, int qlen);
interval  exactMatch( const esa_t *C, const char *query);

interval getInterval( const esa_t *C, const interval ij, char a);
lcp_inter_t getLCPInterval( const esa_t *C, const char *query, size_t qlen);

#endif

