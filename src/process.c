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
 * @param C - The enhanced suffix array of the subject.
 * @param query - The actual query string.
 * @param ql - The length of the query string. Needed for speed reasons.
 */
double dist( esa_t *C, char *query, size_t ql){
	int jumps = 0;
	size_t idx = 0;
	
	
	while( idx < ql ){
		saidx_t l = longestMatch( C, query + idx, ql - idx);
		if( l == 0 ) break;

		jumps++;
		idx += l + 1; // skip the mutation
	}

	// avoid NaN
	if( jumps <= 1){
		jumps = 2;
	}
	return (double)(jumps-1)/(double)ql;
}


double dist_inc( esa_t *C, const char *query, size_t ql){
	size_t jumps = 0; // The jumps found so far
	size_t homol = 0; // Number of homologous nucleotides so far.
	
	size_t projected = 0;
	lcp_inter_t inter;
	
	size_t idx = 0;
	size_t found;
	size_t l;
	
	while( idx < ql ){
		inter = getLCPInterval( C, query + idx, ql - idx );
		l = inter.l;
		if( l == 0){
			break;
		}
		
		if( l >= 100 ){
			fprintf( stderr, "WTF? l:%ld\n", l);
		}
		
		// lets just assume inter.j == inter.i
		found = C->SA[ inter.j];
		
		if( projected == found ){
			// We have homology
			jumps++;
			homol += l + 1;
			idx += l + 1;
			projected = found + l + 1;
		} else {
			idx += l + 1;
			projected = found + l + 1;
		}
	}	
	
	// avoid NaN
	if( jumps <= 1){
		jumps = 2;
	}
	if( homol <= 0){
		homol = 1;
	}
	
	if( jumps >= homol){
		fprintf( stderr, "jumps: %ld, homol: %ld\n", jumps, homol);
		fprintf( stderr, "idx: %ld, ql:%ld, l: %ld\n", idx, ql, l);
		return 0.74999; // 0.75 - epsilon
	}

	return (double)(jumps -1)/(double)homol ;
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

			size_t ql = sequences[j].len;
			
			if( STRATEGY == S_SIMPLE ){
				d = dist( &E, sequences[j].S, ql);
			} else if( STRATEGY == S_INC ){
				d = dist_inc( &E, sequences[j].S, ql);
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

