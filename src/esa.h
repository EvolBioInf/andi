#ifndef _ESA_H_
#define _ESA_H_

#include <divsufsort64.h>
#include <RMQ.hpp>


typedef struct {
	const char *S;
	saidx64_t *SA;
	saidx64_t *ISA;
	saidx64_t *LCP;
	saidx64_t len;
	RMQ *rmq_lcp;
} esa_t;

typedef struct {
	saidx64_t i;
	saidx64_t j;
} interval;

typedef struct {
	saidx64_t l, i,j;
} lcp_inter_t;


int compute_SA( esa_t *c);
int compute_LCP( esa_t *c);
saidx64_t longestMatch( const esa_t *C, const char *query, int qlen);
interval  exactMatch( const esa_t *C, const char *query);

interval getInterval( const esa_t *C, const interval ij, char a);
lcp_inter_t getLCPInterval( const esa_t *C, const char *query, size_t qlen);

#endif

