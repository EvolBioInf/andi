#include <stdlib.h>
#include <RMQ.hpp>
#include <string.h>
#include <stdio.h>
#include <assert.h>
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
		C->SA = (saidx_t*) malloc(C->len * sizeof(saidx_t));
		if( C->SA == NULL){
			return 2;
		}
	}
	
	saidx_t result;
	result = divsufsort((const unsigned char*)C->S, C->SA, C->len);
	
	return result;
}

/**
 * Computes the LCP and ISA given SA and S. This function uses the method from Kasai et al.
 * @param C The enhanced suffix array.
 * @return 0 iff sucessful
 */
int compute_LCP( esa_t *C){
	const char *S = C->S;
	saidx_t *SA  = C->SA;
	saidx_t len  = C->len;
	
	if( !C || S == NULL || SA == NULL || len == 0){
		return 1;
	}
	
	if( C->ISA == NULL){
		C->ISA = (saidx_t*) malloc( len*sizeof(saidx_t));
		if( C->ISA == NULL ){
			return 2;
		}
	}
	if( C->LCP == NULL){
		C->LCP = (saidx_t*) malloc((len+1)*sizeof(saidx_t));
		if( C->LCP == NULL ){
			return 3;
		}
	}
	
	saidx_t *ISA = C->ISA;
	saidx_t *LCP = C->LCP;

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
 * This function implements an alternative way of computing an LCP
 * array for a given suffix array. It uses an intermediate `phi`
 * array, hence the name. It's a bit faster than the other version.
 * @param {esa_t*} C - The enhanced suffix array to compute the LCP from.
 * @returns 0 iff successful
 */
int compute_LCP_PHI( esa_t *C){
	const char *S = C->S;
	saidx_t *SA  = C->SA;
	saidx_t len  = C->len;
	
	// Trivial safety checks
	if( !C || S == NULL || SA == NULL || len == 0){
		return 1;
	}
	
	// Allocate new memory
	if( C->LCP == NULL){
		// The LCP array is one element longer than S.
		C->LCP = (saidx_t*) malloc((len+1)*sizeof(saidx_t));
		if( C->LCP == NULL ){
			return 3;
		}
	}
	saidx_t *LCP = C->LCP;
	
	LCP[0] = -1;
	LCP[len] = -1;
	
	// Allocate temporary arrays
	saidx_t *PHI = (saidx_t *) malloc( len * sizeof(saidx_t));
	saidx_t *PLCP = (saidx_t *) malloc( len * sizeof(saidx_t));
	if( !PHI || !PLCP) return 2;
	
	PHI[SA[0]] = -1;
	ssize_t i, k;
	
	
	for( i=1; i< len; i++){
		PHI[SA[i]] = SA[ i-1];
	}
	
	ssize_t l = 0;
	for( i = 0; i< len ; i++){
		k = PHI[i];
		if( k != -1 ){
			while( S[k+l] == S[i+l]){
				l++;
			}
			PLCP[i] = l;
			l--;
			if( l < 0) l = 0;
		} else {
			PLCP[i] = -1;
		}
	}
	
	// unpermutate the LCP array
	for( i=1; i< len; i++){
		LCP[i] = PLCP[SA[i]];
	}
	
	free(PHI);
	free(PLCP);
	return 0;
}


/**
 * Given the LCP interval for a string `w` this function calculates the 
 * LCP interval for `wa` where `a` is a single character.
 * @param {const esa_t*} C - This is the enhanced suffix array of the subject sequence.
 * @param {lcp_inter_t*} ij - The prefix `w` is given implicitly by the ij LCP interval.
 * @param {char} a - The next character in the query sequence.
 * @returns A reference to the new LCP interval.
 */
lcp_inter_t *getInterval( const esa_t *C, lcp_inter_t *ij, char a){
	saidx_t i = ij->i;
	saidx_t j = ij->j;


	const saidx_t *SA = C->SA;
	const saidx_t *LCP = C->LCP;
	const char *S = C->S;
	RMQ *rmq_lcp = C->rmq_lcp;
	
	if( i == j ){
		if( S[SA[i]] == a){
			ij->i = ij->j = i;
		} else {
			ij->i = ij->j = -1;
		}
		return ij;
	}

	saidx_t l, m;
	
	m = rmq_lcp->query(i+1, j); // m is now any minimum in (i..j]
	l = LCP[m];
	ij->l = l;
	
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
		ij->i = i;
		ij->j = j;
		ij->l = LCP[m];
	} else {
		ij->i = ij->j = -1;
	}
	
	return ij;
}

/* Ohlebusch getInterval Alg 5.2 p.119
 */
lcp_inter_t getLCPInterval( const esa_t *C, const char *query, size_t qlen){
	lcp_inter_t res = {0,0,0};

	if( !C || !query || !C->len || !C->SA || !C->LCP || !C->S || !C->rmq_lcp ){
		res.i = res.j = res.l = -1;
		return res;
	}
	
	saidx_t k = 0, l, i, j, t, p;
	lcp_inter_t ij = { 0, 0, C->len-1};
	saidx_t m = qlen;
	
	saidx_t *SA = C->SA;
	const char *S = (const char *)C->S;
	
	
	do {
		getInterval( C, &ij, query[k]);
		i = ij.i;
		j = ij.j;
		
		if( i == -1 && j == -1 ){
			res.l = k;
			return res;
		}
		
		res.i = ij.i;
		res.j = ij.j;

		l = m;
		if( i < j){
			/* Instead of making another RMQ we can use the LCP interval calculated in getInterval */
			t = ij.l;
			if( t < l ){
				l = t;
			}
		}
		p = SA[i];

		for(;k<l && S[p+k] && query[k];k++){
			if( S[p+k] != query[k] ){
				res.l = k;
				return res;
			}
		}

		k = l;
	} while ( k < m);

	res.l = m;
	return res;
}

/**
 * @returns longest prefix of query found in subject
 */
saidx_t longestMatch( const esa_t *C, const char *query, int qlen){
	return getLCPInterval( C, query, qlen).l;
}

/**
 * exactMatch Ohlebusch Alg 5.2
 * @returns LCP Interval for query or [-1..-1] if not found
 */
interval exactMatch( const esa_t *C, const char *query){
	interval ij;
	lcp_inter_t kl = getLCPInterval( C, query, strlen(query));
	
	ij.i = kl.i;
	ij.j = kl.j;
	
	return ij;
}




















