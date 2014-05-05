#include <stdlib.h>
#include <RMQ.hpp>
#include <string.h>
#include <stdio.h>
#include "esa.h"

/**
 * Computes the SA given a string S.
 * @param C The enhanced suffix array to use. Reads C->S, fills C->SA.
 * @return 0 iff successful
 */
int compute_SA(esa_t *C){
	// assert c.S
	if( !C || C->S == NULL){
		return 1;
	}
	if( C->SA == NULL){
		C->SA = (saidx64_t*) malloc(C->len * sizeof(saidx64_t));
		if( C->SA == NULL){
			return 2;
		}
	}
	
	saidx64_t result;
	result = divsufsort64((const unsigned char*)C->S, C->SA, C->len);
	
	return result;
}

/**
 * Computes the LCP and ISA given SA and S.
 * @param C The enhanced suffix array.
 * @return 0 iff sucessful
 */
int compute_LCP( esa_t *C){
	const char *S = C->S;
	saidx64_t *SA  = C->SA;
	saidx64_t len  = C->len;
	
	if( !C || S == NULL || SA == NULL || len == 0){
		return 1;
	}
	
	if( C->ISA == NULL){
		C->ISA = (saidx64_t*) malloc( len*sizeof(saidx64_t));
		if( C->ISA == NULL ){
			return 2;
		}
	}
	if( C->LCP == NULL){
		C->LCP = (saidx64_t*) malloc((len+1)*sizeof(saidx64_t));
		if( C->LCP == NULL ){
			return 3;
		}
	}
	
	saidx64_t *ISA = C->ISA;
	saidx64_t *LCP = C->LCP;

	LCP[0] = -1;
	LCP[len] = -1;
	
	int i,j,k,l;
	for( i=0; i< len; i++){
		ISA[SA[i]] = i;
	}
	
	l=0;
	for( i=0; i< len; i++){
		j = ISA[i];
		if( j> 0) {
			k = SA[j-1];
			while( S[k+l] == S[i+l] ){
				l++;
			}
			LCP[j] = l;
			l--;
			if (l<0) l = 0;
		}
	}
	return 0;
}


/**
 * Ohlebusch Alg 5.1
 */
interval getInterval( const esa_t *C, interval ij, char a){
	interval ret;
	saidx64_t i = ij.i;
	saidx64_t j = ij.j;
	if( !C || !C->S || !C->SA || !C->LCP || !C->rmq_lcp ){
		ij.i = ij.j = -2;
		return ij;
	}

	saidx64_t *SA = C->SA;
	saidx64_t *LCP = C->LCP;
	const char *S = C->S;
	RMQ *rmq_lcp = C->rmq_lcp;
	
	if( i == j ){
		if( S[SA[i]] == a){
			ret.i = ret.j = i;
		} else {
			ret.i = ret.j = -1;
		}
		return ret;
	}

	saidx64_t l, m;
	
	m = rmq_lcp->query(i+1, j); // m is now any minimum in (i..j]
	l = LCP[m];
	do {
		if( S[ SA[m] + l] <= a ){
			i = m;
		} else {
			j = m - 1;
		}
		if( i == j ){
			break;
		}
		
		m = rmq_lcp->query(i+1, j);
	} while( LCP[m] == l);

	if( S[SA[i] + l] == a){
		ret.i = i;
		ret.j = j;
	} else {
		ret.i = ret.j = -1;
	}
	
	return ret;
}

/**
 * exactMatch Ohlebusch Alg 5.2
 * @returns LCP Interval for query or [-1..-1] if not found
 */
interval exactMatch( const esa_t *C, const char *query){
	interval ij;
	if( !C || !query || !C->len || !C->SA || !C->LCP || !C->S || !C->rmq_lcp ){
		ij.i = ij.j = -2;
		return ij;
	}

	saidx64_t *SA = C->SA;
	saidx64_t *LCP = C->LCP;
	const char *S = (const char *)C->S;
	RMQ *rmq_lcp = C->rmq_lcp;
	
	
	saidx64_t k = 0, l, i, j, t, p;
	ij.i = 0;
	ij.j = C->len-1;
	saidx64_t m = strlen( query);
	do {
		ij = getInterval( C, ij, query[k]);
		i = ij.i;
		j = ij.j;
		if( i == -1 && j == -1 ){
			return ij;
		}
		l = m;
		if( i < j){
			t = LCP[ rmq_lcp->query(i+1, j)];
			if( t < l ){
				l = t;
			}
		}
		p = SA[i] + k ;
		if( strncmp(S+p, query+k, l-k) != 0){
			ij.i = ij.j = -1;
			return ij;
		}
		k = l;
	} while ( k < m);
	
	return ij;
}


/**
 * @returns longest prefix of query found in subject
 */
saidx64_t longestMatch( const esa_t *C, const char *query, int qlen){
	if( !C || !query || !C->len || !C->SA || !C->LCP || !C->S || !C->rmq_lcp ){
		return -1;
	}
	
	saidx64_t k = 0, l, i, j, t, p;
	interval ij = { 0, C->len-1};
	saidx64_t m = qlen;
	
	saidx64_t *SA = C->SA;
	saidx64_t *LCP = C->LCP;
	const char *S = (const char *)C->S;
	RMQ *rmq_lcp = C->rmq_lcp;
	
	do {
		ij = getInterval( C, ij, query[k]);
		i = ij.i;
		j = ij.j;
		if( i == -1 && j == -1 ){
			return k;
		}

		l = m;
		if( i < j){
			t = LCP[ rmq_lcp->query(i+1, j)];
			if( t < l ){
				l = t;
			}
		}
		p = SA[i];
		
		for(;k<l && S[p+k] && query[k];k++){
			if( S[p+k] != query[k] ){
				return k;
			}
		}
		
		k = l;
	} while ( k < m);

	return m;
}

/* Ohlebusch getInterval Alg 5.1 p.118
 */
lcp_inter_t getLCPInterval( const esa_t *C, const char *query, size_t qlen){
	lcp_inter_t res = {0,0,0};

	if( !C || !query || !C->len || !C->SA || !C->LCP || !C->S || !C->rmq_lcp ){
		res.i = res.j = res.l = -1;
		return res;
	}
	
	saidx64_t k = 0, l, i, j, t, p;
	interval ij = { 0, C->len-1};
	saidx64_t m = qlen;
	
	saidx64_t *SA = C->SA;
	saidx64_t *LCP = C->LCP;
	const char *S = (const char *)C->S;
	RMQ *rmq_lcp = C->rmq_lcp;
	
	
	do {
		res.i = ij.i;
		res.j = ij.j;
		ij = getInterval( C, ij, query[k]);
		i = ij.i;
		j = ij.j;
		
		if( i == -1 && j == -1 ){
			res.l = k;
			return res;
		}

		l = m;
		if( i < j){
			t = LCP[ rmq_lcp->query(i+1, j)];
			if( t < l ){
				l = t;
			}
		}
		p = SA[i] + k;

		for(;k<l && S[SA[i]+k] && query[k];k++){
			if( S[SA[i]+k] != query[k] ){
				res.l = k;
				return res;
			}
		}

		k = l;
	} while ( k < m);

	res.l = m;
	return res;
}























