/**
 * @file
 * @brief This file contains various distance methods.
 */
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
double dist( const esa_t *C, char *query, size_t ql){
	size_t jumps = 0;
	size_t idx = 0;
	saidx_t l;
	lcp_inter_t inter;
	
	while( idx < ql ){
		inter = getLCPInterval( C, query + idx, ql - idx);
		//saidx_t l = longestMatch( C, query + idx, ql - idx);
		l = inter.l;
		if( l == 0 ) break;
		
		if( FLAGS & F_EXTRA_VERBOSE ){
			fprintf( stderr, "idx: %ld, l: %d\n", idx, l);
		}
		
		if( FLAGS & F_EXTRA_VERBOSE && l == 423 ){
			fprintf( stderr, "[%d..%d]\n", inter.i, inter.j);
			fprintf( stderr, "%s\n", query + idx );
			fprintf( stderr, "%d\n", C->SA[inter.i] );
			fprintf( stderr, "%s\n", C->S + 992); 
		}

		jumps++;
		idx += l + 1; // skip the mutation
	}

	// avoid NaN
	if( jumps <= 1){
		jumps = 2;
	}
	
	if( FLAGS & F_VERBOSE ){
		fprintf( stderr, "jumps: %ld, homol: %ld\n", jumps, ql);
	}
	
	return (double)(jumps-1)/(double)ql;
}

/**
 * @brief Calculates the log_2 of a given integer.
 */
inline size_t log2( size_t num){
	size_t res = 0;
	while( num >>= 1){
		res++;
	}
	return res;
}

/**
 * @brief Calculates the log_4 of a given integer.
 */
inline size_t log4( size_t num){
	return log2( num) >> 1;
}

double dist_anchor( const esa_t *C, const char *query, size_t query_length){
	size_t snps = 0;
	size_t homo = 0;
	
	lcp_inter_t inter;
	
	size_t last_pos_Q = 0;
	size_t last_pos_S = 0;
	size_t last_length = 0;
	size_t last_was_right_anchor = 0;
	
	size_t this_pos_Q = 0;
	size_t this_pos_S;
	size_t this_length;
	
	size_t num_right_anchors = 0;
	
	// Iterate over the complete query.
	while( this_pos_Q < query_length){
		inter = getLCPInterval( C, query + this_pos_Q, query_length - this_pos_Q);
		this_length = inter.l;
		if( this_length == 0) break;
		
		/* TODO: evaluate the result of different conditions */
		if( inter.i == inter.j  
			&& this_length >= 2 * log4(query_length) )
		{
			// We have reached a new anchor.
			this_pos_S = C->SA[ inter.i];
			
			// Check if this can be a right anchor to the last one.
			if( this_pos_Q - last_pos_Q == this_pos_S - last_pos_S ){
				num_right_anchors++;
			
				// Count the SNPs in between.
				size_t i;
				for( i= 0; i< this_pos_Q - last_pos_Q; i++){
					if( C->S[ last_pos_S + i] != query[ last_pos_Q + i] ){
						snps++;
					}
				}
				homo += this_pos_Q - last_pos_Q;
				last_was_right_anchor = 1;
			} else {
				if( last_was_right_anchor){
					// If the last was a right anchor, but with the current one, we 
					// cannot extend, then add its length.
					homo += last_length;
				}
				
				last_was_right_anchor = 0;
			}
			
			// Cache values for later
			last_pos_Q = this_pos_Q;
			last_pos_S = this_pos_S;
			last_length= this_length;
		}
		
		// Advance
		this_pos_Q += this_length + 1;
	}
	
	// We might miss a few nucleotides if the last anchor was also a right anchor.
	if( last_was_right_anchor ){
		homo += last_length;
	}
	
	// Avoid NaN.
	if ( snps <= 2) snps = 2;
	if ( homo <= 3) homo = 3;
	
	if( FLAGS & F_VERBOSE ){
		fprintf( stderr, "snps: %lu, homo: %lu\n", snps, homo);
		fprintf( stderr, "number of right anchors: %lu\n", num_right_anchors);
	}
	
	return (double)snps/(double)homo;
}

double dist_sophisticated( const esa_t *C, const char *query, size_t ql){
	size_t idx = 0;
	size_t snps = 0;
	size_t homo = 0;
	lcp_inter_t inter;
	saidx_t l;
	ssize_t projected = 0;
	ssize_t i;
	ssize_t found;
	int extendable = 0;
	size_t snpt = 0, homt = 0;
	
	size_t uniqs = 0;
	size_t extensions = 0;
	size_t nonextensions = 0;
	size_t adjacent = 0;
	
	while( idx < ql){
		inter = getLCPInterval( C, query + idx, ql - idx);
		l = inter.l;
		if( l == 0) break;
		
		if( FLAGS & F_EXTRA_VERBOSE ){
			fprintf(stderr, "idx: %3ld, inter.i: %4d, found: %4d, %d\n", idx, inter.i, C->SA[inter.i], extendable);
		}
		
		// unique
		if( inter.i == inter.j){
			uniqs++;
			
			found = C->SA[inter.i];
			
			if( extendable && found >= projected && found - projected < l ){
				snps += snpt + 1;
				homo += homt + l + 1;
				adjacent++;
			} else {
				snps += 2;
				homo += l + 2;
			}
			
			idx += l + 1;
			
			projected = found + l + 1;
			extendable = 1; 
		} else if( extendable ){
			// neighbour?
			found = -1;
			
			if( inter.j - inter.i < 50 ){
				for( i= inter.i; i <= inter.j; i++){
					if( C->SA[i] < projected ){
						continue;
					}
				
					if( found == -1 || C->SA[i] < (size_t)found){
						found = (ssize_t) C->SA[i];
					}
				}
			} else {
				found = inter.i;
			}
		
			if( found >= projected && found - projected < l ){
				// we have homology
				snps += snpt + 1;
				homo += homt + l + 1;
				projected = found + l + 1;
				extensions++;
			} else {
				extendable = 0;
				nonextensions++;
			}
			
			idx += l + 1;
		} else {
			nonextensions++;
			idx += l + 1;
		}
		
		snpt = 0;
		homt = 0;
				
		if( extendable ){
			int k = 0; // Number of allowed mismatches
			int LOOKAHEAD = 0;
			for( i= 0; i< LOOKAHEAD && k; i++){
				if( C->S[ projected + i] != query[ idx + i] ){
					idx += i + 1;
					projected += i + 1;
					snpt++;
					homt = i + 1;
					k--;
				}
			}
		}
		
	}
	
	if( FLAGS & F_VERBOSE ){
		fprintf( stderr, "snps: %ld, homo: %ld\n", snps, homo);
		fprintf( stderr, "unique matches: %ld, adjacent: %ld, extensions: %ld, nonextensions: %ld\n", 
			uniqs, adjacent, extensions, nonextensions);
	}
	
	return (double)snps/(double)homo;
}

/**
 * @brief Computes the distance matrix.
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

	#pragma omp parallel for num_threads( THREADS)
	for(i=0;i<n;i++){
		esa_t E = {NULL,NULL,NULL,NULL,0,NULL};
		
		// initialize the enhanced suffix array
		if( FLAGS & F_SINGLE ){
			E.S = (const char*) sequences[i].S;
			E.len = sequences[i].len;
		} else {
			E.S = (const char*) sequences[i].RS;
			E.len = sequences[i].RSlen;
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
			} else if( STRATEGY == S_SOPHISTICATED){
				d = dist_sophisticated( &E, sequences[j].S, ql);
			} else if( STRATEGY == S_ANCHOR ){
				d = dist_anchor( &E, sequences[j].S, ql);
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
 * @brief Prints the distance matrix.
 * @param sequences An array of pointers to the sequences.
 * @param n The number of sequences.
 */
void printDistMatrix( seq_t* sequences, int n){
	int i, j;

	// initialise the sequences
	#pragma omp parallel for num_threads( THREADS)
	for( i=0;i<n;i++){
		if( !sequences[i].S){
			fprintf(stderr, "missing sequence %d\n", i);
			exit(1);
		}
		sequences[i].len = strlen( sequences[i].S);
		if( !(FLAGS & F_SINGLE) ){
			sequences[i].RS = catcomp( sequences[i].S, sequences[i].len);
			sequences[i].RSlen = 2 * sequences[i].len + 1;
		}
	}
	
	double *D = distMatrix( sequences, n);
	
	printf("%d\n", n);
	for( i=0;i<n;i++){
		printf("%8s", sequences[i].name);
		for( j=0;j<n;j++){
			printf(" %1.4lf", D(i,j));
		}
		printf("\n");
	}
	
	free(D);
}


