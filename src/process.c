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
		
		if( FLAGS & F_VERBOSE ){
			fprintf( stderr, "idx: %ld, l: %d\n", idx, l);
		}
		
		if( l == 423 ){
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


double dist_check( const esa_t *C, const char *query, size_t ql){
	size_t idx = 0;
	size_t snps = 0;
	size_t homo = 0;
	lcp_inter_t inter;
	saidx_t l;
	
	size_t i;
	size_t found;
	
	while( idx < ql){
		inter = getLCPInterval( C, query + idx, ql - idx);
		l = inter.l;
		if( l == 0) break;
		
		snps += 1;
		idx += l + 1;
		homo += l + 1;
		
		// lets just assume inter.i == inter.j
		found = C->SA[inter.i] + l + 1;
		for( i=0; i< 2; i++){
			if( C->S[found + i] != query[idx + i]){
				idx += i + 1;
				snps++;
				homo += i + 1;
				break;
			}
		}
	}
	
	if( FLAGS & F_VERBOSE ){
		fprintf( stderr, "snps: %ld, homo: %ld\n", snps, homo);
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

double dist_inc( const esa_t *C, const char *query, size_t ql){
	size_t jumps = 0; // The jumps found so far
	size_t homol = 0; // Number of homologous nucleotides so far.
	
	ssize_t projected = 0;
	lcp_inter_t inter;
	
	size_t idx = 0;
	ssize_t found;
	size_t l;
	
	while( idx < ql ){
		inter = getLCPInterval( C, query + idx, ql - idx );
		l = inter.l;
		if( l == 0){
			break;
		}
		
		if( inter.j - inter.i < 10 ){
			found = -1;
			saidx_t k = inter.i;
			
			for( ; found == -1 && k <= inter.j; k++){
				if( C->SA[k] >= projected ){
					found = C->SA[k];
					break;
				}
			}
			
			for( ; k <= inter.j; k++){
				if( C->SA[k] >= projected && C->SA[k] < found ){
					found = C->SA[k];
				}
			}
			
			if( found == -1 ){
				found = C->SA[ inter.i];
			}
			
		} else {
			found = C->SA[ inter.j];
		}
		
		if( found >= projected && (size_t)(found - projected) < l ){
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
		return 0.74999; // 0.75 - epsilon
	}

	if( FLAGS & F_VERBOSE ){
		fprintf( stderr, "jumps: %ld, homol: %ld\n", jumps, homol);
	}

	return (double)(jumps -1)/(double)homol ;
}


#define WINDOW 5

double dist_window( const esa_t *C, const char *query, size_t ql){
	size_t jumps = 0; // number of jumps so far
	size_t homol = 0; // number of homologous nucleotides so far
	
	lcp_inter_t inter;
	size_t l; // inter.l
	
	ssize_t projected[WINDOW] = {0};
	ssize_t found[WINDOW] = {0};
	size_t  dist[WINDOW] = {0};
	ssize_t candidate;
	size_t k, p;
	size_t count = 0;
	
	size_t idx = 0;
	while( idx < ql ){
		inter = getLCPInterval( C, query + idx, ql - idx );
		l = inter.l;
		if( l == 0){
			break;
		}
		
		if( inter.j - inter.i < 5 ){ // the difference between 2, 5 and 50 is negligible
			for( p=0; p< WINDOW; p++){
				found[p] = -1;
			}
			
			for( k= inter.i; k<= inter.j; k++){
				candidate = C->SA[ k];
				for( p=0; p< WINDOW; p++){
					if( candidate < projected[p]){
						continue;
					}
					
					if( found[p] == -1 || candidate < found[p]){
						found[p] = candidate;
					}
				}
			}
			
			for( p=0; p< WINDOW; p++){
				if( found[p] == -1 ){
					found[p] = inter.i;
				}
			}
		} else {
			for( p=0; p< WINDOW; p++){
				found[p] = inter.i;
			}
		}
		
		for( p=0; p< WINDOW; p++){
			if( found[p] >= projected[p] && (size_t)(found[p] - projected[p]) < l ){
				// homology
				jumps += p + 1;
				homol += dist[p] + l + 1;
				if( p > 0 ) count ++;
				break;
			}
		}
		
		for( p= WINDOW; p-- > 1; ){
			dist[p] = dist[p-1] + l + 1;
			projected[p] = projected[p-1] + l + 1;
		}
		
		dist[0] = 0;
		projected[0] = found[0] + l + 1;
		
		idx += l + 1;
	}
	
	
	// avoid NaN
	if( jumps <= 1){
		jumps = 2;
	}
	if( homol <= 0){
		homol = 1;
	}
	
	if( jumps >= homol){
		return 0.74999; // 0.75 - epsilon
	}

	if( FLAGS & F_VERBOSE ){
		fprintf( stderr, "jumps: %ld, homol: %ld\n", jumps, homol);
		fprintf( stderr, "jumps with history >= 1: %ld\n", count);
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
			} else if( STRATEGY == S_INC ){
				d = dist_inc( &E, sequences[j].S, ql);
			} else if( STRATEGY == S_WINDOW){
				d = dist_window( &E, sequences[j].S, ql);
			} else if( STRATEGY == S_CHECK){
				d = dist_check( &E, sequences[j].S, ql);
			} else if( STRATEGY == S_SOPHISTICATED){
				d = dist_sophisticated( &E, sequences[j].S, ql);
			}
			
			if( !(FLAGS & F_RAW)){
				/*  Our shustring method might miss a mutation or two. Hence we need to correct  
					the distance using math. See Haubold, Pfaffelhuber et al. (2009) */
				if( STRATEGY == S_SIMPLE || STRATEGY == S_INC 
					|| STRATEGY == S_WINDOW || STRATEGY == S_CHECK ){
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
		if( !(FLAGS & F_SINGLE) ){
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


