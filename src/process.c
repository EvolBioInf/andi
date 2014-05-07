#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "esa.h"
#include "divergence.h"
#include "global.h"
#include "process.h"

#include <RMQ_1_n.hpp>
#include <RMQ_succinct.hpp>


#define D( X, Y) (D[ (X)*n + (Y) ])

/**
 * The dist function computes the distance between a subject and a query
 * @return Number of jumps which is related to the number of mutations
 * @param C The enhanced suffix array of the subject.
 * @param query The actual query string.
 * @param ql The length of the query string. Needed for speed reasons.
 */
int dist( esa_t *C, char *query, int ql){
	int jumps = 0;
	saidx_t idx = 0;
	
	
	while( idx < ql ){
		saidx_t l = longestMatch( C, query + idx, ql - idx);
		if( l == 0 ) break;

		jumps++;
		idx += l + 1; // skip the mutation
	}

	return jumps;
}


int dist_inc( esa_t *C, char *query, int ql){
	int jumps = 0;
	int extral = 0;
	int extrar = 0;
	saidx_t idx = 0;
	saidx_t expected = 0;
	lcp_inter_t inter = {0,0,0};
	saidx_t lastpos, thispos;
	saidx_t found;
	RMQ *rmq_SA = new RMQ_succinct( C->SA, C->len);
	saidx_t *SA = C->SA;
	saidx_t k;
	signed long sldist;
	
	double dist = 0;
	lastpos = 0;
	
	while( idx < ql ){
		inter = getLCPInterval( C, (char*) (query+idx), ql-idx);
		if( inter.l == 0 ) break;
		
		
		if( inter.j - inter.i < 4 ){
			found = -1;
			k = inter.i;
			
			for(; k<= inter.j && found == -1; k++){
				if( SA[k] >= expected ){
					found = k;
					break;
				}
			}
			
			for( ; k <= inter.j; k++){
				if( SA[k] >= expected && SA[k] < SA[found]){
					found = k;
				}
			}
			if( found == -1 ){
				found = inter.i;
			}
		} else {
			found = rmq_SA->query(inter.i, inter.j);
		}
		thispos = SA[ found];
		
		// Assume thispos > lastpos. Otherwise we skipped a mutation (either here or there).
		sldist = (long)thispos - (long)lastpos;
		dist += fabs((double)sldist);
		
		if( FLAGS & F_EXTRA_VERBOSE){
			fprintf(stderr, "%d\t%d\t%d\t%2d-[%d..%d]\n", idx, expected, thispos, inter.l, inter.i, inter.j);
		}
		
		if( -sldist > (signed long)(2 * (inter.l + 1)) ){
			// we propably missed another mutation
			extral++;
			expected += inter.l + 1;
		} else if( sldist > (signed long)(2 * (inter.l + 1)) ){
			extrar++; // dito
			expected += inter.l + 1;
		} else {
			// no extra mutation
			expected = thispos + inter.l + 1;
		}
		lastpos = thispos;
		
		jumps++;
		idx += inter.l + 1; // skip the mutation
	}
	
	if( FLAGS & F_VERBOSE ){
		printf("jumps: %d, extral: %d, extrar: %d\n", jumps, extral, extrar);
		printf("var dist: %lf\n", dist/(double)jumps);
	}

	return jumps + extral + extrar ;
}


/**
 * The distMatrix populates the D matrix with computed distances. It allocates D and
 * filles it with useful values, but the caller has to free it!
 * @return The distance matrix
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
double *distMatrix( seq_t* sequences, int n){
	double *D = (double*)malloc( n * n * sizeof(double));
	assert(D);
	
	double d;
	
	int i;

	#pragma omp parallel for num_threads( CORES)
	for(i=0;i<n;i++){
		esa_t E = {NULL,NULL,NULL,NULL,0,NULL};
		
		// initialize the enhanced suffix array
		if( FLAGS & F_DOUBLE ){
			E.S = (const char*) sequences[i].RS;
			E.len = sequences[i].RSlen;
		} else {
			E.S = (const char*) sequences[i].S;
			E.len = sequences[i].len;
		}
		
		int result;

		result = compute_SA( &E);
		assert( result == 0); // zero errors
		result = compute_LCP_PHI( &E);
		assert( result == 0);
	
		E.rmq_lcp = new RMQ_succinct(E.LCP, E.len);
		
		// now compare every other sequence to i
		int j;
		for(j=0; j<n; j++){
			if( j == i) {
				D(i,j) = 0.0;
				continue;
			}
			
			if( FLAGS & F_VERBOSE ){
				#pragma omp critical
				{
					fprintf( stderr, "comparing %d and %d\n", i, j);
				}
			}

			int ql = sequences[j].len;
			
			if( STRATEGY == S_SIMPLE ){
				result = dist( &E, sequences[j].S, ql);
			} else if( STRATEGY == S_INC ){
				result = dist_inc( &E, sequences[j].S, ql);
			}
			
			if( result <= 2) result = 3; // avoid NaN
			
			if( FLAGS & F_VERBOSE ){
				printf("i: %d, j: %d, jumps: %d, length: %d\n", i, j, result, ql);
			}
			
			
			if( FLAGS & F_DOUBLE ){
				d = (double)(result - 2)/(double)ql;
			} else {
				d = (double)(result - 1)/(double)ql;
			}
			
			if( !(FLAGS & F_RAW)){
				/*  Our shustring method might miss a mutation or two. Hence we need to correct  
					the distance using math. See Haubold, Pfaffelhuber et al. (2009) */
				if( STRATEGY == S_SIMPLE ){
					d = divergence( 1.0/d, E.len, 0.5, 0.5);
				}
				d = -0.75 * log(1.0- (4.0 / 3.0) * d ); // jukes cantor
			}
			D(i,j) = d;
		}
		
		delete E.rmq_lcp;
		free( E.SA);
		free( E.ISA);
		free( E.LCP);
	}
	
	return D;
}

/**
 * Prints the distance matrix.
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
void printDistMatrix( seq_t* sequences, int n){
	int i, j;

	// initialise the sequences
	#pragma omp parallel for num_threads( CORES)
	for( i=0;i<n;i++){
		if( !sequences[i].S){
			fprintf(stderr, "missing sequence %d\n", i);
			exit(1);
		}
		sequences[i].len = strlen( sequences[i].S);
		if( FLAGS & F_DOUBLE ){
			sequences[i].RS = catcomp( sequences[i].S, sequences[i].len);
			sequences[i].RSlen = 2 * sequences[i].len + 1;
		}
	}
	
	double *D = distMatrix( sequences, n);
	
	for( i=0;i<n;i++){
		printf("%8s", sequences[i].name);
		for( j=0;j<n;j++){
			printf(" %1.4lf", D(i,j));
		}
		printf("\n");
	}
	
	free(D);
}



/**
 * Compute the reverse complement.
 * @param str The master string.
 * @param len The length of the master string
 * @return The reverse complement. The caller has to free it!
 */
char *revcomp( const char *str, size_t len){
	char *rev = (char*) malloc( len + 1);
	if( !str || !rev) return NULL;
	
	char *r = rev;
	const char *s = str + len-1;
	
	rev[len] = '\0';
	
	char c, d;
	
	while( len --> 0 ){
		c = *s--;
		
		switch( c){
			case 'A': d = 'T'; break;
			case 'T': d = 'A'; break;
			case 'G': d = 'C'; break;
			case 'C': d = 'G'; break;
			default: continue;
		}
		
		*r++ = d;
	}
	
	return rev;
}

/**
 * This function concatenates the reverse complement to a given master string.
 * @param s - The master string.
 * @param len - Its length.
 * @return The new concatenated string.
 */
char *catcomp( char *s , size_t len){
	if( !s) return NULL;
	
	char *rev = revcomp( s, len);
	
	rev = (char*) realloc( rev, 2 * len + 2);
	if( !rev) return NULL;
	
	rev[len] = '#';
	
	memcpy( rev+len+1, s, len+1);
	
	return rev;
}

