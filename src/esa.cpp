/**
 * @file
 * @brief ESA functions
 *
 * This file contains various functions that operate on an enhanced suffix
 * array. Some of these are taken straight from the book of Ohlebusch (2013);
 * Others are modified for improved performance. For future reference here are
 * some of the implemented changes.
 *
 * The name `getLCPInterval` is kind of missleading. It and `getInterval` both
 * find LCP-intervals but once for prefixes and for characters respectively.
 * A critical speed component for both functions is the number of RMQs done.
 * To reduce the number of calls, the `m` property of `lcp_inter_t` is caching
 * the correct value. This reduces the number of RMQ calls per input base
 * to an average of 1.7. This is less than the expected number of calls for
 * a binary search over for elements (log_2 (4) = 2).
 *
 * An additional RMQ-saving strategy was introduced in commit 7ba1d363189… and
 * optimized in 716958a01…. Since for every new match the very same RMQs are
 * called a cache is introduced: A "hash" over the first `CACHE_LENGTH` characters
 * of a new query determines an index into the cache. This cache is filled at
 * the initilization time of the ESA. Hence for multiple comparisions this is
 * done only once. The filling itself is done via a depth first search over the
 * imaginary suffix tree. This strategy saves about 20% of time when comparing
 * only two sequences and provides an additional speedup a factor of 3 to 4 when
 * applied to big data sets.
 */
#include <stdlib.h>
#include <string.h>
#include "esa.h"

static void esa_init_cache_dfs( esa_t *C, char *str, size_t pos, const lcp_inter_t *in);
void esa_init_cache_fill( esa_t *C, char *str, size_t pos, const lcp_inter_t *in);

lcp_inter_t *getNoRMQInterval( const esa_t *C, lcp_inter_t *ij, char a);
lcp_inter_t getNoRMQLCPInterval( const esa_t *C, const char *query, size_t qlen);
lcp_inter_t getNoRMQLCPIntervalFrom( const esa_t *C, const char *query, size_t qlen, saidx_t k, lcp_inter_t ij);

static int esa_init_SA( esa_t *c);
static int esa_init_LCP( esa_t *c);
int esa_init_CLD( esa_t *C);

/** @brief The prefix length up to which RMQs are cached. */
const size_t CACHE_LENGTH = 8;

/** @brief Map a code to the character. */
char code2char( ssize_t code){
	switch( code & 0x3){
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
	}
	return '\0';
}

/** @brief Map a character to a two bit code. */
ssize_t char2code( const char c){
	ssize_t result = -1;
	switch( c){
		case 'A' : result = 0; break;
		case 'C' : result = 1; break;
		case 'G' : result = 2; break;
		case 'T' : result = 3; break;
	}
	return result;
}

/** @brief Fills the RMQ cache.
 *
 * When looking up LCP intervals of matches we constantly do the same RMQs on
 * the LCP array. Since a RMQ takes about 45 cycles we can speed up things by
 * caching the intervals for some prefix length (`CACHE_LENGTH`).
 *
 * @param C - The ESA.
 * @returns 0 iff successful
 */
int esa_init_cache( esa_t *C){
	lcp_inter_t* rmq_cache = (lcp_inter_t*) malloc((1 << (2*CACHE_LENGTH)) * sizeof(lcp_inter_t) );

	if( !rmq_cache){
		return 1;
	}

	C->rmq_cache = rmq_cache;

	char str[CACHE_LENGTH+1];
	str[CACHE_LENGTH] = '\0';
	lcp_inter_t ij = { 0, 0, C->len-1, 0};
	//ij.m = C->rmq_lcp->query(1,C->len-1);
	esa_init_cache_dfs( C, str, 0, &ij);

	return 0;
}

/** @brief Fills the RMQ cache — one char at a time.
 *
 * Recursively traverse the virtual sufix tree created by the SA, LCP and RMQs,
 * using a depth first search. by doing this cache the current LCP interval and
 * improve the performance of getting the next value for the cache.
 */
void esa_init_cache_dfs( esa_t *C, char *str, size_t pos, const lcp_inter_t *in){
	// we are not yet done, but the current strings do not exist in the subject.
	if( pos < CACHE_LENGTH && in->i == -1 && in->j == -1){
		esa_init_cache_fill(C,str,pos,in);
		return;
	}

	if( pos >= CACHE_LENGTH){
		esa_init_cache_fill(C,str,pos,in);
		return;
	}

	lcp_inter_t ij;

	for( int code = 0; code < 4; ++code){
		str[pos] = code2char(code);
		ij = *in;
		getNoRMQInterval(C, &ij, str[pos]);

		// fail early
		if( ij.i == -1 && ij.j == -1){
			esa_init_cache_fill(C, str, pos + 1, &ij);
			continue;
		}

		if ( ij.l > (ssize_t)(pos + 1)){

			if( (size_t)ij.l < CACHE_LENGTH){
				// fill with dummy value
				esa_init_cache_fill(C, str, pos+1, in);

				// fast forward
				auto k = pos + 1;
				for(;k < (size_t)ij.l; k++){
					str[k] = C->S[C->SA[ij.i]+k];
				}
				esa_init_cache_dfs(C, str, k, &ij);
				continue;
			}

			esa_init_cache_fill(C, str, pos+1, in);
			continue;
		}

		esa_init_cache_dfs(C,str,pos+1,&ij);
	}
}

void esa_init_cache_fill( esa_t *C, char *str, size_t pos, const lcp_inter_t *in){
	if( pos < CACHE_LENGTH){
		for( int code = 0; code < 4; ++code){
			str[pos] = code2char(code);
			esa_init_cache_fill( C, str, pos + 1, in);
		}
	} else {
		size_t code = 0;
		for( size_t i = 0; i < CACHE_LENGTH; ++i ){
			code <<= 2;
			code |= char2code(str[i]);
		}

		C->rmq_cache[code] = *in;
	}
}

/** @brief Initializes an ESA.
 *
 * This function initializes an ESA with respect to the provided sequence.
 * @param C - The ESA to initialize.
 * @param S - The sequence
 * @returns 0 iff successful
 */
int esa_init( esa_t *C, seq_t *S){
	C->S = NULL;
	C->SA = NULL;
	C->ISA = NULL;
	C->LCP = NULL;
	C->len = 0;
	C->CLD = NULL;

	int result;

	C->S = (const char*) S->RS;
	C->len = S->RSlen;

	if( C->S == NULL ) return 1;

	result = esa_init_SA(C);
	if(result) return result;

	result = esa_init_LCP(C);
	if(result) return result;

	result = esa_init_CLD(C);
	if(result) return result;

	result = esa_init_cache(C);
	if(result) return result;

	return 0;
}

/** @brief Free the private data of an ESA. */
void esa_free( esa_t *C){
	free( C->SA);
	free( C->ISA);
	free( C->LCP);
	free( C->rmq_cache);
	free( C->CLD);
}

/**
 * Computes the SA given a string S. To do so it uses libdivsufsort.
 * @param C The enhanced suffix array to use. Reads C->S, fills C->SA.
 * @returns 0 iff successful
 */
int esa_init_SA(esa_t *C){
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

#define R(CLD, i) ((CLD)[(i)])
#define L(CLD, i) ((CLD)[(i)-1])

int esa_init_CLD( esa_t *C){
	if( !C || !C->LCP){
		return 1;
	}
	saidx_t* CLD =  C->CLD = (saidx_t*) malloc((C->len+1) * sizeof(saidx_t));
	if( !C->CLD) {
		return 2;
	}

	saidx_t *LCP = C->LCP;

	typedef struct pair_s {
		saidx_t idx, lcp;
	} pair_t;

	pair_t *stack = (pair_t*) malloc((C->len+1) * sizeof(pair_t));
	pair_t *top = stack;
	pair_t last;

	R(CLD,0) = C->len + 1;

	top->idx = 0;
	top->lcp = -1;
	top++;
	ssize_t k = 1;
	for(; k < C->len + 1; k++){
		while( LCP[k] < (top-1)->lcp){
			last = *(--top);

			while( (top-1)->lcp == last.lcp){
				R(CLD,(top-1)->idx) = last.idx;
				last = *--top;
			}

			if( LCP[k] < (top-1)->lcp){
				R(CLD, (top-1)->idx) = last.idx;
			} else {
				L(CLD, k) = last.idx;
			}
		}
		top->idx = k;
		top->lcp = LCP[k];
		top++;
	}

	free( stack);
	return 0;
}

/**
 * This function implements an alternative way of computing an LCP
 * array for a given suffix array. It uses an intermediate `phi`
 * array, hence the name. It's a bit faster than the other version.
 * @param C The enhanced suffix array to compute the LCP from.
 * @returns 0 iff successful
 */
int esa_init_LCP( esa_t *C){
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
	if( !PHI || !PLCP){
		free(PHI);
		free(PLCP);
		return 2;
	}
	
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

lcp_inter_t *getNoRMQInterval( const esa_t *C, lcp_inter_t *ij, char a){
	int i = ij->i;
	int j = ij->j;
	int m;
	

	const saidx_t *SA = C->SA;
	const saidx_t *LCP = C->LCP;
	const char *S = C->S;
	const saidx_t *CLD = C->CLD;

	// check for singleton or empty interval
	if( i == j ){
		if( S[SA[i]] == a){
			ij->i = ij->j = i;
		} else {
			ij->i = ij->j = -1;
		}
		return ij;
	}

	if( LCP[i] <= LCP[j+1]){
		m = L(CLD, j+1);
	} else {
		m = R(CLD, i);
	}

	ssize_t k = ij->i;
	int l = LCP[m];


	do {
		if( S[ SA[k] + l] == a ){
			ij->i = k;
			ij->j = m-1;
			ij->m = LCP[ij->i] <= LCP[ij->j+1] ? L(CLD, m) : R(CLD,k) ;
			ij->l = LCP[ij->m];
			//ij->l = LCP[C->rmq_lcp->query(k+1, m-1)]; //LCP[]
			return ij;
		}
		k = m;
		if( k == j){
			break;
		}
		m = R(CLD,m);
	} while ( /*m != "bottom" && */ LCP[m] == l);

	// final sanity check
	if( S[SA[k] + l] == a){
		ij->i = k;
		ij->j = j;
		/* Also return the length of the LCP interval including `a` and
		 * possibly even more characters. Note: l + 1 <= LCP[m] */
		ij->l = LCP[m];
		ij->m = m;
	} else {
		ij->i = ij->j = -1;
	}

	return ij;
}

lcp_inter_t getNoRMQLCPIntervalFrom( const esa_t *C, const char *query, size_t qlen, saidx_t k, lcp_inter_t ij){

	if( ij.i == -1 && ij.j == -1){
		return ij;
	}

	// fail early on singleton intervals.
	if( ij.i == ij.j){

		// try to extend the match. See line 499 below.
		saidx_t p = C->SA[ij.i];
		size_t k = ij.l;
		const char *S = (const char *)C->S;

		for(k = 0 ; k< qlen && S[p+k]; k++ ){
			if( S[p+k] != query[k]){
				ij.l = k;
				return ij;
			}
		}

		ij.l = k;
		return ij;
	}

	saidx_t l, i, j, p;
	saidx_t m = qlen;

	lcp_inter_t res = ij;
	
	saidx_t *SA = C->SA;
	const char *S = (const char *)C->S;
	
	// Loop over the query until a mismatch is found
	do {
		getNoRMQInterval( C, &ij, query[k]);
		i = ij.i;
		j = ij.j;
		
		// If our match cannot be extended further, return.
		if( i == -1 && j == -1 ){
			res.l = k;
			return res;
		}
		
		res.i = ij.i;
		res.j = ij.j;

		l = m;
		if( i < j){
			/* Instead of making another RMQ we can use the LCP interval calculated
			 * in getInterval */
			if( ij.l < l ){
				l = ij.l;
			}
		}
		
		// Extend the match
		for( p = SA[i]; k < l && S[p+k] && query[k]; k++){
			if( S[p+k] != query[k] ){
				res.l = k;
				return res;
			}
		}
		
		// TODO: Verify if this is the best solution
		// You shall not pass the null byte.
		if( k < l && (!S[p+k] || !query[k])){
			res.l = k;
			return res;
		}

		k = l;
	} while ( k < m);

	res.l = m;
	return res;
}


lcp_inter_t getNoRMQLCPInterval( const esa_t *C, const char *query, size_t qlen){
	lcp_inter_t res = {0,0,0,0};

	// sanity checks
	if( !C || !query || !C->len || !C->SA || !C->LCP || !C->S || !C->CLD ){
		res.i = res.j = res.l = -1;
		return res;
	}

	lcp_inter_t ij = { 0, 0, C->len-1, 0};

	return getNoRMQLCPIntervalFrom(C, query, qlen, 0, ij);
}


/** @brief Compute the LCP interval of a query. For a certain prefix length of the
 * query its LCP interval is retrieved from a cache. Hence this is faster than the
 * naive `getInterval`.
 *
 * @param C - The enhanced suffix array for the subject.
 * @param query - The query sequence.
 * @param qlen - The length of the query. Should correspond to `strlen(query)`.
 * @returns The LCP interval for the longest prefix.
 */
lcp_inter_t getNoRMQCachedLCPInterval( const esa_t *C, const char *query, size_t qlen){
	if( qlen <= CACHE_LENGTH) return getNoRMQLCPInterval( C, query, qlen);

	ssize_t offset = 0;
	for( size_t i = 0; i< CACHE_LENGTH; i++){
		offset <<= 2;
		offset |= char2code(query[i]);
	}

	if( offset < 0){
		return getNoRMQLCPInterval( C, query, qlen);
	}

	lcp_inter_t ij = C->rmq_cache[offset];

	if( ij.j == -1 && ij.j == -1){
		return getNoRMQLCPInterval( C, query, qlen);
	}

	return getNoRMQLCPIntervalFrom(C, query, qlen, ij.l, ij);
}
